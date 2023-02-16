clc
clear all
close all

mtot = 1580;                                        %[kg]
mu = 160;                                           %[kg]
ks1 = 27*10^3; ks2 = 22.5*10^3; k1p = 4*274*10^3;   %[N/m]
c = 2664;   %value of 1DOF
ms = mtot-mu;
k = 2*(ks1+ks2);

c_opt = sqrt(ms*k/2*(k1p+2*k)/k1p);

M = [ms 0; 
    0 mu];

K = [k -k; 
    -k k+k1p];

C = [c -c; 
    -c c];

H = [0 k1p]';

%% Modes evaluation
[vv,ww] = eig(K,M); % produces a diagonal matrix ww
    % of generalized eigenvalues and a full matrix vv
    % whose columns are the corresponding eigenvectors
    % so that kk*vv = mm*vv*ww --> (kk-mm*ww)*vv = 0
om = diag(ww).^0.5 % natural frequencies
    % eigenvector normalisation in order to have the
    % first element of phi = 1
f_nat = om/2/pi
phi_1 = vv(1:2,1)./(vv(1,1));
phi_2 = vv(1:2,2)./(vv(1,2));
phi = [phi_1 , phi_2]

%% Compute TF with receptance method
OM = [0:0.1:180].';
C_values = [c_opt c_opt/2 c_opt/3 2*c_opt];
C_names = {'c_opt' 'c_opt/2' 'c_opt/3' '2*c_opt'};

for count2=1:length(C_values)
    c = C_values(count2);
    C = [c -c;-c c];
    for count1=1:size(OM,1)
        w = OM(count1,1);
        Kdyn = (K + i*w.*C - w^2.*M);
        Rec = inv(Kdyn);
        X12(count1,1) = Rec(1,2)*H(2,1); % tf: zs/h
        X22(count1,1) = Rec(2,2)*H(2,1); % tf: zu/h
        Fz(count1,1) = k1p*(1-X22(count1,1)); % tf: Fz/h
    end
    figure(1); 
    subplot(211); hold all;
    plot(OM/(2*pi),abs(X12),'-.','Displayname',C_names{count2});
    ylabel('X_{12} = z_s/h)'); legend('show'); grid on;
    
    subplot(212); hold all; 
    plot(OM/(2*pi),abs(X22),'-','Displayname',C_names{count2});
    ylabel('X_{22} = z_u/h'); xlabel('frequency [Hz]')
    legend('show'); grid on;
    
    figure(2); 
    subplot(211); hold all;
    plot(OM/(2*pi),OM.^2.*abs(X12),':','Displayname',C_names{count2});
    ylabel('zpp_s/h'); legend('show'); grid on;
    
    subplot(212); hold all;
    plot(OM/(2*pi),abs(Fz)./k1p,'--','Displayname',C_names{count2});
    ylabel('F_z/(k_p h)'); xlabel('frequency [Hz]')
    legend('show'); grid on;
end