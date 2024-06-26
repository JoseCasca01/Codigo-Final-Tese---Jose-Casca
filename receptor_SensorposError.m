function [val,Rerror] = receptor_SensorposError(fieldx,fieldy,BS,R,group,f,c,N,lambda,erro,mode)

    %The receiver will wait a time interval equal to t for all signals
    
    PT = 1;
    GT = 1;
    GR = 1;
    
    PR = zeros(N,1);
    for i = 1:N
        %Friis formula for free-space propagation    
        PR(i) = PT * GT * GR * lambda^2/(16*pi^2*R(i)^2);
    end
    
    %The sensor furthest from the receiver will be the reference sensor with phase 0
    phaseCorrect = (max(R)-R).*2*pi/lambda;
    traveling_time = R/c;

    PropagationLoss = (4*pi*R/lambda).^2;
    grouperror = zeros(size(group,1),size(group,2),size(group,3));
    for i = 1:size(group,1)
        for j = 1:size(group,2)
                %randn de média = 0 e variancia = erro
                grouperror(i,j) = group(i,j) + (erro*randn(1,1));
        end
    end
    
    %To remove z component error from sensors -> activate line below
    %grouperror(:,:,3)=zeros(size(grouperror,1),size(grouperror,2),1);

    Rerror = distance(grouperror,BS);    
    t=linspace(0,max(Rerror)*3/c,10000); % 1 período

    phase = (max(Rerror)-Rerror).*2*pi/lambda;
    
    signals = zeros(N,length(t));
    signalsCorrect = zeros(N,length(t));

    for i = 1:N
        signals(i,:) = cos(2*pi*f*t-phase(i));
        signalsCorrect(i,:) = cos(2*pi*f*t-phaseCorrect(i));
    end
    
    [~, adjustsignals_Loss, receivedsignal, receivedsignal_Loss] = adjsig(N,t,signals,traveling_time,PropagationLoss);
    [~, ~, ~, idealsignal_Loss] = adjsig(N,t,signalsCorrect,traveling_time,PropagationLoss);

    aux = sum(adjustsignals_Loss);
    val = max(aux(t>=1.1*traveling_time(R==max(R))));

    if mode == 1
        close all;

        figure(1);
        %%%%% 2D%%%%%
        plot(BS(1),BS(2),'X'),hold on;
        plot(group(:,1),group(:,2),'O')
        plot(grouperror(:,1),grouperror(:,2),'diamond');
        title('Field');
        ylabel('yfield (m)');
        xlabel('xfield (m)');
        legend('Drone', 'Node Position', 'Node Position with error');
        axis([0, fieldy, 0, fieldx]);
        %%%%%3D%%%%%
        figure(8);
        plot3(BS(1),BS(2),BS(3),'X'),hold on;
        plot3(group(:,1),group(:,2),group(:,3),'O')
        plot3(grouperror(:,1),grouperror(:,2),grouperror(:,3),'diamond');
        title('Field');
        ylabel('yfield (m)');
        xlabel('xfield (m)');
        legend('Drone', 'Node Position', 'Node Position with error');
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
        plot(t,receivedsignal);
        title('Received Signal');
        xlabel('t(s)');
        ylabel('Amplitude');
        axis([0 max(t)*1.01 min(receivedsignal)-1 max(receivedsignal)+1]);
        
        % figure(5);
        % plot(t,adjustsignals_Loss);
        % title('Signals Received with Losses With Position Error');
        % xlabel('t(s)');
        % ylabel('Amplitude');
        
        figure(6);
        plot(t,receivedsignal_Loss);
        title('Received Signal with Losses');
        xlabel('t(s)');
        ylabel('Amplitude');
        axis([0 max(t)*1.01 min(receivedsignal_Loss)*1.01 max(receivedsignal_Loss)*1.01]);

        figure(7);
        plot(t,idealsignal_Loss-receivedsignal_Loss);
        title('Received Signals Sum Differences');
        xlabel('t(s)');
        ylabel('Amplitude');
        
        %Signals plot without time adjust
        %figure(10);
        %plot(t,signals(:,:));

        disp(val)
    end

end