function R = distance(transmitters,receiver)
    % transmitters(N,2) -> N elements with 2 dimensions
    % receiver(1,1) -> receiver position
    % R -> (N,1) matrix that saves the distance between transmitter i and
    % receiver
    N = size(transmitters,1);
    R = zeros(N,1);
    for i = 1:N
        R(i) = sqrt((transmitters(i,1)-receiver(1))^2+(transmitters(i,2)-receiver(2))^2+(transmitters(i,3)-receiver(3))^2);
    end
end