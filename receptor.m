function [val] = receptor(R,f,c,N,lambda,mode)
    %AF = 0;
    %O primeiro sensor a transmitir será o mais distante ao recetor. Este será
    %o sensor de referência com fase 0
    phase = (max(R)-R).*2*pi/lambda;
    traveling_time = R/c;
    
    %O receptor vai esperar um intervalo de tempo igual a t pelos sinais todos
    tmax = max(R)*3/c;
    t=linspace(0,tmax,10000); % 1 período
    %tstem = 0:1/f:tmax;

    %signals = ones(N,length(t));
    signals = zeros(N,length(t));
    for i = 1:N
        signals(i,:) = cos(2*pi*f*t-phase(i));
        %AF = AF + exp(1j*2*pi*(-R(i))/lambda);
    end
    
    PT = 1;
    GT = 1;
    GR = 1;
    
    PR = zeros(N,1);
    for i = 1:N
        %Friis formula for free-space propagation    
        PR(i) = PT * GT * GR * lambda^2/(16*pi^2*R(i)^2);
    end
    
    PropagationLoss = (4*pi*R/lambda).^2;

    [adjustsignals, adjustsignals_Loss, receivedsignal,receivedsignal_Loss] = adjsig(N,t,signals,traveling_time,PropagationLoss);

    %Signals plot without time adjust
    %figure(10);
    %plot(t,signals(:,:));
    if mode == 1
        figure(2);
        plot(1:N,PropagationLoss,'-X');
        title('Propagation Loss per sensor')
        ylabel('P_R');
        xlabel('Sensor number');
        
        figure(3);
        plot(t,adjustsignals(:,:));
        title('Received Signals');
        xlabel('t(s)');
        ylabel('Amplitude');
        
        figure(4);
        plot(t,receivedsignal);
        title('Received Signal');
        xlabel('t(s)');
        ylabel('Amplitude');
        axis([0 max(t)*1.01 min(receivedsignal)-1 max(receivedsignal)+1]);
        
        figure(5);
        plot(t,adjustsignals_Loss);
        title('Received Signals with Losses');
        xlabel('t(s)');
        ylabel('Amplitude');
        
        figure(6);
        plot(t,receivedsignal_Loss);
        title('Received Signal with Losses');
        xlabel('t(s)');
        ylabel('Amplitude');
        axis([0 max(t)*1.01 min(receivedsignal_Loss)*1.01 max(receivedsignal_Loss)*1.01]);
        %figure(7);
        %plot(tstem,receivedsignal(t == tstem));
    end

    val = max(sum(adjustsignals_Loss));

end