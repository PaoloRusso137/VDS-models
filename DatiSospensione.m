clc
clear vars
close all

c = 2664; 
k = 35495; %N/m
ms = (1580-160)/4;
mu = 160/4;
m = ms+mu;
ex_freq = 10; % Road excitation frequency [rad/s]
h_amp = 0.01; % Ampltoude road elevation [m]  
kp = 274e03; %tyre stiffness
Ts = 4;

% nfile = 'Simscape_1DOF';
% open(nfile)
% sim(nfile)
% 
% 
% plot(Velocity)
% xlabel('Time {s}')
% ylabel('Velocity [m/s]')
% 
% plot(Positioning)
% xlabel('Time {s}')
% ylabel('Positioning [m]')