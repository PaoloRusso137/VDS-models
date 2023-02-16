%% SIMULATION SET UP
clc
close all

dati = 'Dati_Necessari';
open(dati)
run(dati)

% mtot = 1580;                                        %[kg]
% mu = 160;                                           %[kg]
% ks1 = 27*10^3; ks2 = 22.5*10^3; k1p = 4*274*10^3;   %[N/m]
% c = 2664;   %value of 1DOF
% ms = mtot-mu;
% k = 2*(ks1+ks2);
% 
% c_opt = sqrt(ms*k/2*(k1p+2*k)/k1p);

choice1 = menu('Values of C','c = 2664','c optimal')

switch choice1
    case 2
        c = c_opt;
%     case 3 
%         c = c_cr;
end

M = [ms 0; 
    0 mu];

K = [k -k; 
    -k k+k1p];

C = [c -c; 
    -c c];

H = [0 k1p]';

% sim setup
z0 = 0;
zd0 = 0;

% Step attributes
lambda = 2;
h = 5/100; %5 cm = 5/100 m
t0 = 2;
slope = h/0.01; %slope for the step signal
v = 0;

choice2 = menu('Road profile h','Step signal h0 = 5cm','Sinewave','Chirp')

switch choice2
            case 1
                Ts=10;
            case 2
                Ts=10;
                choice3 = menu('Select the velocity','V = 10 km/h','V = 50 km/h','V = 130 km/h')
                switch choice3
                    case 1
                        v = 10/3.6;
                    case 2
                        v = 50/3.6;
                    case 3
                        v = 130/3.6;
                end
            case 3
                Ts=10;
end

if choice2==1
    Ts=10;
end

if choice2==3
    Ts=100;
end

w_in = 0;           %[rad/s]
w_end = 180;        %[rad/s]
f_in=0.1;     %[Hz]
f_end=w_end/2/pi;   %[Hz]
h0 = 2/100;
f_sampling = 1*10^3; T_sampling = 1/f_sampling;
decimation = 1;

nfile = 'sim_2DOF_ANIM_TF';
open(nfile)
sim(nfile)

%% Animation & Outputs for sin & step

prompt = "Should animation start? Y/N [Y]: ";
txt = input(prompt,"s");
if txt == 'Y'
    decim = 100;
    ti = 0;
    tf = Ts;
    anim_2DOF(v,road_profile,z,time,decim,ti,tf)
    return
end

if choice2==1 || choice2==2
    figure(1)
    plot(zdd(:,1)) %zdd from simulink is in m!!!
    xlabel('time [s]'); xlim([0 1]);
    ylabel('Sprung mass Acceleration [m/s^2]')
    legend('Sprung mass','Unsprung mass'); grid on; hold on

    figure(2)
    plot(time,zs_NL,'LineWidth',1.5) %zdd from simulink is in m!!!
    xlabel('time [s]'); xlim([0 Ts]);
    ylabel('Sprung mass Displacement [m/s^2]')
    grid on
    hold on

    figure(2)
    plot(time,zu_NL,'LineWidth',1.5) %zdd from simulink is in m!!!
    xlabel('time [s]'); xlim([0 Ts]);
    ylabel('Unprung mass Displacement [m/s^2]')
    legend('Sprung mass','Unsprung mass'); grid on; hold on

    return
end

%% TF estimate
[G_zdds_h,F] = tfestimate(road_profile,zdd.data(:,1),[],[],[],f_sampling); %zs is the first column of zdd

relative_disp = road_profile-zu;
G_Fz_kph = tfestimate(k1p*road_profile,k1p*relative_disp);

f = linspace(f_in,f_end,length(G_zdds_h));

figure(1) % SPRUNG ACCELERATION
plot(F,abs(G_zdds_h)) 
% semilogx(F,20*log10(abs(G_zdds_h)))
grid on
xlabel('Frequency [Hz]')
ylabel('Sprung mass acceleration TF')
xlim([f_in f_end])
hold on

figure(2) % VERTICAL FORCE TF
plot(F,abs(G_Fz_kph)) 
% semilogx(F,20*log10(abs(G_Fz_kph)))
grid on
xlabel('Frequency [Hz]')
ylabel('Vertical force TF')
xlim([f_in f_end])
legend('Optimal damping')
hold on

%% Non linear damping 

v_NL_damp = [-1000 -520 -390 -260 -130 -50 0 50 130 260 390 520 1000]/1000;
F_NL_damp = 4*[-672 -493 -450 -400 -351 -199 0 301 1068 1260 1386 1514 2130];

nfile = 'sim_NONLINEAR_2DOF';
open(nfile)
sim(nfile)

[G_zdds_h_NL,F] = tfestimate(road_profile_NL,zdd_NL.data(:,1),[],[],[],f_sampling); %zs is the first column of zdd

relative_disp_NL = road_profile_NL-zu_NL;
G_Fz_kph_NL = tfestimate(k1p*road_profile_NL,k1p*relative_disp_NL);

f = linspace(f_in,f_end,length(G_zdds_h_NL));

% TRANSFER FUNCTIONS
figure(1) % SPRUNG ACCELERATION
plot(F,abs(G_zdds_h_NL)) 
% semilogx(F,20*log10(abs(G_zdds_h)))
grid on
xlabel('Frequency [Hz]')
ylabel('Sprung mass acceleration TF')
xlim([f_in f_end])
legend('Optimal damping','Non linear damping')

figure(2) % VERTICAL FORCE TF
plot(F,abs(G_Fz_kph_NL)) 
% semilogx(F,20*log10(abs(G_Fz_kph)))
grid on
xlabel('Frequency [Hz]')
ylabel('Vertical force TF')
xlim([f_in f_end])
legend('Optimal damping','Non linear damping')