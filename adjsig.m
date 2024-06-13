function [newsignals,newsignals_Loss, receivedsignal, receivedsignal_Loss] = adjsig(N,t,signals,trav_time,PropagationLoss)
    
    newsignals=zeros(N,length(t));
    newsignals_Loss = zeros(N,length(t));

    for i = 1:N
        newsignals(i,t>trav_time(i)) = signals(i,1:end-length(t(t<=trav_time(i))));
        newsignals_Loss(i,:) = newsignals(i,:) / PropagationLoss(i);
    end
    
    receivedsignal = sum(newsignals);
    receivedsignal_Loss = sum(newsignals_Loss);
end