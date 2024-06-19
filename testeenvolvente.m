close all;

f=1e3;
w0=2*pi;

t=[0.25 0.75];
ak=sin(w0*t);

t1=0:0.001:1;
ak1 = sin(2*pi*t1);

figure;
plot(t,ak),hold on;
plot(t1,ak1);

%%
figure;
plot(t1,ak1);