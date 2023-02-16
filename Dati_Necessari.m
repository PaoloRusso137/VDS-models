warning('off','all')

%% MATLAB DATA
mtot = 1580;                                        %[kg]
mu = 160;                                           %[kg]
ks1 = 27*10^3; ks2 = 22.5*10^3; k1p = 4*274*10^3;   %[N/m]
c = 2664;   %value of 1DOF
ms = mtot-mu;
k = 2*(ks1+ks2);

c_opt = sqrt(ms*k/2*(k1p+2*k)/k1p);

M = [ms 0; 0 mu];

K = [k -k; -k k+k1p];

C = [c -c; -c c];

H = [0 k1p]';

%% SIMULINK DATA
% Initial conditions
z0 = 0;
zd0 = 0;

% Chirp attributes
w_in = 0;           %[rad/s]
w_end = 180;        %[rad/s]
f_in=0.1;     %[Hz]
f_end=w_end/2/pi;   %[Hz]
h0 = 2/100;
f_sampling = 1*10^3; T_sampling = 1/f_sampling;
decimation = 1;

% Step attributes
lambda = 2;
h = 5/100; %5 cm = 5/100 m
t0 = 2;
slope = h/0.01; %slope for the step signal
v = 0;