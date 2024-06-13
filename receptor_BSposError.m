function [val] = receptor_BSposError(fieldx,fieldy,BS,R,group,f,c,N,lambda,erro,mode)

    %O receptor vai esperar um intervalo de tempo igual a t pelos sinais todos
    t=linspace(0,max(R)*3/c,10000); % 1 período
    
    PT = 1;
    GT = 1;
    GR = 1;
    
    PR = zeros(N,1);
    for i = 1:N
        %Friis formula for free-space propagation    
        PR(i) = PT * GT * GR * lambda^2/(16*pi^2*R(i)^2);
    end
    
    %O primeiro sensor a transmitir será o mais distante ao recetor. Este será
    %o sensor de referência com fase 0
    phaseCorrect = (max(R)-R).*2*pi/lambda;
    traveling_time = R/c;

    PropagationLoss = (4*pi*R/lambda).^2;
    BSerror = BS + (erro*randn(1,1));
    
    Rerror = distance(group,BSerror);    
    
    phase = (max(Rerror)-Rerror).*2*pi/lambda;
    
    signals = zeros(N,length(t));
    signalsCorrect = zeros(N,length(t));

    for i = 1:N
        signals(i,:) = cos(2*pi*f*t-phase(i));
        signalsCorrect(i,:) = cos(2*pi*f*t-phaseCorrect(i));
    end

    [adjustsignals,adjustsignals_Loss] = adjsig(N,t,signals,traveling_time,PropagationLoss);
    [~, adjustsignals_LossCorrect] = adjsig(N,t,signalsCorrect,traveling_time,PropagationLoss);
    %Valor espera < 0 visto que o sinal com erro na posição terá -x Potência do que o sinal ideal 

    aux = sum(adjustsignals_Loss);
    val=max(aux(t>=1.1*traveling_time(R==max(R))));
    
    if mode == 1
        close all;

        figure(1);
        plot(group(:,1),group(:,2),'O'), hold on;
        plot(BS(1),BS(2),'X');
        plot(BSerror(1),BSerror(2),'diamond');
        title('Field');
        ylabel('yfield (m)');
        xlabel('xfield (m)');
        axis([0, fieldy, 0, fieldx]);

        figure(2);
        plot(1:N,PropagationLoss,'-X');
        title('Propagation Loss per sensor')
        ylabel('P_R');
        xlabel('Sensor number');
    
        % figure(3);
        % plot(t,adjustsignals(:,:));
        % title('Signals Received With Position Error');
        % xlabel('t(s)');
        % ylabel('Amplitude');
        
        figure(4);
        plot(t,sum(adjustsignals));
        title('Signals Received Sum');
        xlabel('t(s)');
        ylabel('Amplitude');
        
        % figure(5);
        % plot(t,adjustsignals_Loss);
        % title('Signals Received with Losses With Position Error');
        % xlabel('t(s)');
        % ylabel('Amplitude');
        
        figure(6);
        plot(t,sum(adjustsignals_Loss));
        title('Signals Received Sum with Losses');
        xlabel('t(s)');
        ylabel('Amplitude');
    
        figure(7);
        plot(t,sum(adjustsignals_LossCorrect)-sum(adjustsignals_Loss));
        title('Received Signals Sum Differences');
        xlabel('t(s)');
        ylabel('Amplitude');
        
        %Signals plot without time adjust
        %figure(10);
        %plot(t,signals(:,:));

        disp(val)
    end
end