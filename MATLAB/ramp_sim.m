clear all;
close all;
dw = 5e4/5e6;
t = 0:dw:100;
N = length(t);
err = zeros(1, N);
p = zeros(1, N);
kp = -1.25;
ki = -1.35;
integrator = 0;
p(1) = pi;
for i = 1:N-1
    err(i) = p(i)-t(i);
    integrator = integrator + ki * err(i);
    p(i+1) = p(i)+kp*err(i)+integrator;
end

plot(t, err)