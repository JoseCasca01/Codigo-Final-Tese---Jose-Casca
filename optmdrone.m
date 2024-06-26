function [valsignal, BSoptm, Roptm] = optmdrone(fieldx,fieldy,group,f,c,N,lambda,divx,divy)
    valsignal = zeros(divx,divy);
    spacex = fieldx/divx;
    spacey = fieldy/divy;

    BSx = zeros(1,divx);
    BSy = zeros(1,divy);
    
    for i = 0:(divx-1)
        BSx(i+1) = (i*spacex+(i+1)*spacex)/2;
    end
    for i = 0:(divy-1)
        BSy(i+1) = (i*spacey+(i+1)*spacey)/2;
    end
    
    maxsignal = zeros(1,10000);

    for i = 1:divx
        for j = 1:divy

            BS = [BSx(i) BSy(j) 20*lambda];

            R = distance(group,BS);

            phase = (max(R)-R).*2*pi/lambda;
            traveling_time = R/c;
    
            %The receiver will wait a time interval equal to t for all signals
            tmax = max(R)*3/c;
            t=linspace(0,tmax,10000); % 1 período

            %signals = ones(N,length(t));
            signals = zeros(N,length(t));
            for m = 1:N
                signals(m,:) = cos(2*pi*f*t-phase(m));
                %AF = AF + exp(1j*2*pi*(-R(m))/lambda);
            end
    
            PT = 1;
            GT = 1;
            GR = 1;
    
            PR = zeros(N,1);
            for m = 1:N
                %Friis formula for free-space propagation    
                PR(m) = PT * GT * GR * lambda^2/(16*pi^2*R(m)^2);
            end
    
        PropagationLoss = (4*pi*R/lambda).^2;

        [~, ~, ~,receivedsignal_Loss] = adjsig(N,t,signals,traveling_time,PropagationLoss);
        
        valsignal(divy-j+1,i) = max(abs(receivedsignal_Loss));

        if max(receivedsignal_Loss) > max(maxsignal)
            maxsignal = receivedsignal_Loss;
            BSoptm = BS; 
            Roptm = R;
        end

        % figure(1);
        % plot(BS(1),BS(2),'bX');
        % hold on;
        % title('Drone positions');
        % ylabel('yfield (m)');
        % xlabel('xfield (m)');
        % axis([0, fieldy*1.01, 0, fieldx*1.01]);
        % legend('Drone');

        % figure(4);
        % plot(t,receivedsignal),hold on;
        % title('Received Signal');
        % xlabel('t(s)');
        % ylabel('Amplitude');
        % axis([0 max(t)*1.01 min(receivedsignal)-1 max(receivedsignal)+1]);

        end
    end

        figure(2);
        plot(group(:,1),group(:,2),'rO');
        title('Field');
        ylabel('yfield (m)');
        xlabel('xfield (m)');
        axis([0, fieldy*1.01, 0, fieldx*1.01]);
        legend('Sensors');
end
