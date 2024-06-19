close all;
clear;
clc;

f = 3e9; % 3 GHz
c = 3e8;
lambda = c/f;
flag = true;

while flag
    N = input('Number of nodes: ');
    if mod(N,1)==0
        flag = false; 
    end
end

%Define a rectangular field
fieldx = 10*lambda; %0 to fieldx
fieldy = 10*lambda; %0 to fieldy

%Base station coordinates
BS = [fieldx/2,fieldy/2,20*lambda]; %center with z = 5
if BS(3) < 1
    BS(3) = 1;
elseif BS(3)>10
    BS(3) = 10;
end
%Define the number of sensors per group
sensorsPos = zeros(N,3); %rows - sensors; collums - coordenates x and y

%The groups need to have sensors close to each other
%Groups dimension 1.5 by 1.5
%Define groups limits 
%sensorsPos -> 0.5 to 2 on X and 0.5 to 2 on Y

%sensorsPos(1:end,1) = 0.5 + 1.5 * rand(size(sensorsPos,1),1);
%sensorsPos(1:end,2) = 0.5 + 1.5 * rand(size(sensorsPos,1),1);

%Random sensors in all the field
sensorsPos(1:end,1) = fieldx * rand(size(sensorsPos,1),1);
sensorsPos(1:end,2) = fieldy * rand(size(sensorsPos,1),1);

figure(1);
plot3(BS(1),BS(2),BS(3),'rX'), hold on;
plot3(sensorsPos(:,1),sensorsPos(:,2),sensorsPos(:,3),'bO');
title(['Field with ', num2str(N), ' nodes']);
ylabel('yfield (m)');
xlabel('xfield (m)');
legend('Base Station', 'Nodes');
axis([0, fieldy, 0, fieldx]);

%% Receção sem erro na posição
R = distance(sensorsPos,BS);
val = receptor(R,f,c,N,lambda,1);

%% Receção com erro gaussiano na posição dos sensores
%Quando não são apresentadas figuras é apresentada a média da diferença de
%potência da posição ideal para a posição com erro
figures = 1; %0 para não ver as figures / 1 para ver figuras
R = distance(sensorsPos,BS);
if figures ~= 0
    receptor_SensorposError(fieldx,fieldy,BS,R,sensorsPos,f,c,N,lambda,0.2*lambda,figures);
else
    val = zeros(1,1000);
    for i = 1:length(val)
        val(i)=receptor_SensorposError(fieldx,fieldy,BS,R,sensorsPos,f,c,N,lambda,0.2*lambda,figures);
    end
    disp(mean(val));
end

%% Receção com erro gaussiano na posição do recetor
%Quando não são apresentadas figuras é apresentada a média da diferença de
%potência da posição ideal para a posição com erro
figures = 1; %0 para não ver as figures / 1 para ver figuras
R = distance(sensorsPos,BS);

if figures ~= 0
    receptor_BSposError(fieldx,fieldy,BS,R,sensorsPos,f,c,N,lambda,0.2*lambda,figures);
else
    val = zeros(1,1000);
    for i = 1:length(val)
        val(i) = receptor_BSposError(fieldx,fieldy,BS,R,sensorsPos,f,c,N,lambda,0.2*lambda,figures);
    end
    disp(mean(val));
end

%% Signal Amplitude Ideal vs Sensors Erro
k=1;

R = distance(sensorsPos,BS);

val = receptor(R,f,c,N,lambda,0);

rounds = 300;
variance = 0.05:0.025:0.5;
Median_y =  zeros(1,length(variance));
Median_x = zeros(1,length(variance));
valtests = zeros(rounds,length(variance));

for m = 1:1
    for i = 1:length(variance)
        for j = 1:rounds*m
            valtests(j,i) = receptor_SensorposError(fieldx,fieldy,BS,R,sensorsPos,f,c,N,lambda,variance(i)*lambda,0);
        end
    end


    valnormalized = valtests/val;

    figure;
    hd = boxplot(valnormalized,variance);
    hold on;
    xlabel('Variance (wavelength)');
    ylabel('Normalized received signal power');
    title('Received signal normalized by ideal signal');
    for i=6:7:7*length(variance)-1
        Median_x(k) = sum(get(hd(i),'XData'))/2;
        Median_y(k) = sum(get(hd(i),'YData'))/2;
        k=k+1;
    end
    plot(Median_x,Median_y);
    legend('Average');
end

%% Signal Amplitude Ideal vs BS Erro
clear variance
clear valnormalized

k=1;

R = distance(sensorsPos,BS);

val = receptor(R,f,c,N,lambda,0);

rounds = 1000;

%variance1(1,:) = 0.025:0.05:3;
%variance2(1,:) = 3.025:0.05:6;
%variance3(1,:) = 0.025:0.05:0.75;
%variance4(1,:) = 2:0.2:3;
variance = 3.025:0.05:6;
Median_y =  zeros(1,length(variance));
Median_x = zeros(1,length(variance));
valtests = zeros(rounds,length(variance));
  
for m = 1:1
    for i = 1:length(variance)
        for j = 1:rounds
            [valtests(j,i),~] = receptor_BSposError(fieldx,fieldy,BS,R,sensorsPos,f,c,N,lambda,variance(i)*lambda,0);
        end
    end


    valnormalized = valtests/val;

    figure;
    hd = boxplot(valnormalized,variance);
    hold on;
    xlabel('Variance (wavelength)');
    ylabel('Normalized received signal power');
    title('Received signal normalized by ideal signal');
    for i=6:7:7*length(variance)-1
        Median_x(k) = sum(get(hd(i),'XData'))/2;
        Median_y(k) = sum(get(hd(i),'YData'))/2;
        k=k+1;
    end
    plot(Median_x,Median_y);
    legend('Average');
end

%% Array Factor Ideal vs Sensors Erro
clear variance
clear valnormalized
    
k = 1;

R = distance(sensorsPos,BS);

AFIdeal = N; %Com a posição ideal o AF é igual ao número de sensores

rounds = 300;
variance = 0.025:0.025:0.5;
Median_y =  zeros(1,length(variance));
Median_x = zeros(1,length(variance));

for m = 1:1
    AFtest = zeros(rounds,length(variance));
    for i = 1:length(variance)
        for j = 1:rounds*m
            [~, grouperror] = receptor_SensorposError(fieldx,fieldy,BS,R,sensorsPos,f,c,N,lambda,variance(i)*lambda,0);
            Rerror = distance(grouperror,BS);
            AFtest(j,i) = sum(exp(1j*2*pi/lambda*(R-Rerror)));
        end
    end
    AFnormalized = abs(AFtest/AFIdeal);

    figure;
    hd = boxplot(AFnormalized,variance);
    hold on;
    xlabel('Variance (wavelength)');
    ylabel('Normalized Array Factor');
    title('Array Factor normalized by ideal signal')
    for i=6:7:7*length(variance)-1
        Median_x(k) = sum(get(hd(i),'XData'))/2;
        Median_y(k) = sum(get(hd(i),'YData'))/2;
        k=k+1;
    end
    plot(Median_x,Median_y);
    legend('Average');
end

%% Array Factor Ideal vs BS Erro
k = 1;

R = distance(sensorsPos,BS);

AFIdeal = N; %Com a posição ideal o AF é igual ao número de sensores

rounds = 1000;
%variance1(1,:) = 0.025:0.05:3;
%variance2(1,:) = 3.025:0.05:6;
%variance3(1,:) = 0.025:0.05:0.75;
%variance4(1,:) = 2:0.2:3;
variance = 3.025:0.05:6;
Median_y =  zeros(1,length(variance));
Median_x = zeros(1,length(variance));

for m = 1:1
    AFtest = zeros(rounds,length(variance));
    for i = 1:length(variance)
        for j = 1:rounds*m
            [~, Rerror] = receptor_BSposError(fieldx,fieldy,BS,R,sensorsPos,f,c,N,lambda,variance(i)*lambda,0);
            AFtest(j,i) = sum(exp(1j*2*pi/lambda*(R-Rerror)));
        end
    end
    AFnormalized = abs(AFtest/AFIdeal);

    figure;
    hd = boxplot(AFnormalized,variance);
    hold on;
    xlabel('Variance (wavelength)');
    ylabel('Normalized Array Factor');
    title('Array Factor normalized by ideal Array Factor')
    for i=6:7:7*length(variance)-1
        Median_x(k) = sum(get(hd(i),'XData'))/2;
        Median_y(k) = sum(get(hd(i),'YData'))/2;
        k=k+1;
    end
    plot(Median_x,Median_y);
    legend('Average');
end

%% Ponto ótimo do drone
close all;

[val,BSoptm,R] = optmdrone(fieldx,fieldy,sensorsPos,f,c,N,lambda,25,25);

%x = 0:fieldx/25:fieldx;
%y = fieldy:-fieldy/25:0;

figure;
imagesc(val/max(max(val)))
colorMap = jet(256);
colormap(colorMap);
colorbar;

%%
k = 1;


AFIdeal = N; %Com a posição ideal o AF é igual ao número de sensores

rounds = 1000;
%variance1(1,:) = 0.025:0.05:3;
%variance2(1,:) = 3.025:0.05:6;
%variance3(1,:) = 0.025:0.05:0.75;
%variance4(1,:) = 2:0.2:3;
variance = 0.025:0.05:3;
Median_y =  zeros(1,length(variance));
Median_x = zeros(1,length(variance));

for m = 1:1
    AFtest = zeros(rounds,length(variance));
    for i = 1:length(variance)
        for j = 1:rounds*m
            [~, Rerror] = receptor_BSposError(fieldx,fieldy,BSoptm,R,sensorsPos,f,c,N,lambda,variance(i)*lambda,0);
            AFtest(j,i) = sum(exp(1j*2*pi/lambda*(R-Rerror)));
        end
    end
    AFnormalized = abs(AFtest/AFIdeal);

    figure;
    hd = boxplot(AFnormalized,variance);
    hold on;
    xlabel('Variance (wavelength)');
    ylabel('Normalized Array Factor');
    title('Array Factor normalized by ideal Array Factor')
    for i=6:7:7*length(variance)-1
        Median_x(k) = sum(get(hd(i),'XData'))/2;
        Median_y(k) = sum(get(hd(i),'YData'))/2;
        k=k+1;
    end
    plot(Median_x,Median_y);
    legend('Average');
end
