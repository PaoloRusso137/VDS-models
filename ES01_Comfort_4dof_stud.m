clc; close all; clc; warning('off','all')
% 4 DOF Vehicle model for Comfort Analysis

%% Vehicle data 
% vehicle data 
clear; close all; clc

%------------ masses and inertias ----------------------
% sprung mass [kg]
m =  1790;
% unsprung mass - front [kg]
m1 =  90;
% unsprung mass - rear [kg] 
m2 =  70;
% mass moment of inertia I_yy [kg m^2]
Jy = 3.293e3;

%--------------- geometry --------------------------------
% wheelbase [m]
L = 2.885;
% front wheelbase [m]
a = 1.45;
% rear wheelbase [m]
b = L-a;
% % centre of mass height [m]
% h_g = 0.49;

%----------------- suspensions ----------------------------
% front axle stiffness [N/m]
ks1 = 2*29.1e3;
% rear axle stiffness [N/m]
ks2 = 2*31.5e3;
% front axle damping [Ns/m]
cs1 = 2*5.7e3;
% rear axle damping [Ns/m]
cs2 = 2*4.5e3;

%------------------- tyres (225/50 R 17) --------------------
% N.B.: axle stiffness (2*tyre stiffness)
% tyre stiffness [N/m]
k1 = 2*2.26954e5;
k2 = k1;
% tyre damping coefficient [Ns/m]
c1 = 2*50;
c2 = c1;
% % wheel moment of inertia [kg m^2] (2 wheels)
% Jns = 2*0.8;

% optimum damping 
cott = sqrt(m*(ks1+ks2)/2*(k1+k2+2*(ks1+ks2))/(k1+k2));
cs1 = cott/(cs1+cs2)*cs1 /1.2;
cs2 = cott/(cs1+cs2)*cs2 /1.2;

%% System Matrices
% Mass
mm = diag([m,Jy,m1,m2])
% Damping
cc = [cs1+cs2   -cs1*a+cs2*b         -cs1    -cs2
    -cs1*a+cs2*b   cs1*a^2+cs2*b^2   cs1*a   -cs2*b
    -cs1               cs1*a          c1+cs1     0
    -cs2              -cs2*b             0       c2+cs2]
% Stiffness
kk = [ks1+ks2   -ks1*a+ks2*b         -ks1    -ks2
    -ks1*a+ks2*b   ks1*a^2+ks2*b^2   ks1*a    -ks2*b
    -ks1               ks1*a          k1+ks1     0
    -ks2              -ks2*b             0       k2+ks2]

%% Modal analysis of the undamped system
% eigenvalues and eigenvectors 
disp('eigenvetors and eigenvalues {ZG,theta,Z1,Z2}')
[vv,ww] = eig(kk,mm)
% pause

% natural frequencies
om = sqrt(diag(ww));
f_n = om/2/pi

% eigenvectors normalisation (first element of the vectors set equal to 1)
psi_1 = vv(:,1)./(vv(1,1));   % 1st eigenvector
psi_2 = vv(:,2)./(vv(1,2));   % 2nd eigenvector
psi_3 = vv(:,3)./(vv(1,3));   % 3rd eigenvector
psi_4 = vv(:,4)./(vv(1,4));   % 4th eigenvector

% modal matrix
disp('Modal Matrix')
psi = [psi_1 , psi_2, psi_3 , psi_4]

% position of the nodes related to each mode shape
disp('nodes: distance from G')
x1_0 = psi_1(1)/psi_1(2)   %[m] node position of the first mode
x2_0 = psi_2(1)/psi_2(2)   %[m] node position of the second mode
x3_0 = psi_3(1)/psi_3(2)   %[m] node position of the 3rd mode
x4_0 = psi_4(1)/psi_4(2)   %[m] node position of the 4th mode

%% Plot mode shape

%% Mode 1
scale_factor = 2e-1;

h_fig_MS = figure('Name','Modal shapes','units','normalized','outerposition',[0 0.5 0.5 0.5]);
x = -b:0.01:a;     % array of the coordinate along beam length
y1 = scale_factor*(psi_1(1) - x*psi_1(2)); % vertical coordinates of the beam (1st mode)

% plot mode shapes
plot(x,y1,'-k','linewidth',3); hold on % 1st mode
plot([x1_0], [0],'xr','linewidth',2,'markersize',12) % nodes

% extension up to the node of mode 1
x_node = x1_0:0.01:-b;
y1_node = scale_factor*(psi_1(1) - x_node*psi_1(2));
plot(x_node,y1_node,'--b','linewidth',1); hold on

% Plot horizontal line
xlim1 = xlim; ylim1 = ylim;
line(xlim1.',[0; 0],...
    'linewidth',1,...
    'color',[0,0,0],'linestyle',':')
text([x(1),0, x(end)], [y1(1),scale_factor*psi_1(1),y1(end)],{'B','G','A'}, 'HorizontalAlignment','center',...
    'BackgroundColor',[1 1 1],'FontSize',12)
% Plot vertical line
axis equal
title(['1st mode shape (f_1=',num2str(round(f_n(1)*100)/100),' Hz)']); xlabel('x [m]'); ylabel('y'); 


% Plot sprung mass
hold on;
R = 0.1;
offsetZ_m_ns = -0.4;

h3 = rectangle('Position',[-b-R,offsetZ_m_ns-R,2*R,2*R],'FaceColor','none','Curvature',[0,0],'linestyle','--','linewidth',0.5);
h4 = rectangle('Position',[-b-R,offsetZ_m_ns+scale_factor*psi_1(4)-R,2*R,2*R],'FaceColor','b','Curvature',[0.8,0.4]);
h5 = rectangle('Position',[a-1*R,offsetZ_m_ns-R,2*R,2*R],'FaceColor','none','Curvature',[0,0],'linestyle','--','linewidth',0.5);
h6 = rectangle('Position',[a-1*R,offsetZ_m_ns+scale_factor*psi_1(3)-R,2*R,2*R],'FaceColor','b','Curvature',[0.8,0.4]);
grid on
% wheels
R = 0.3;
h3 = rectangle('Position',[-b-R,1*offsetZ_m_ns-R,2*R,2*R],'FaceColor','none','Curvature',[1,1]);
h5 = rectangle('Position',[a-1*R,1*offsetZ_m_ns-R,2*R,2*R],'FaceColor','none','Curvature',[1,1]);

% road
area(xlim,offsetZ_m_ns-R*[1 1]-0.3,offsetZ_m_ns-R,'linewidth',1,'FaceColor',[128 128 128]./256); hold on;


%% Mode 2
scale_factor = 7e-2;

h_fig_MS2 = figure('Name','Modal shapes','units','normalized','outerposition',[0.5 0.5 0.5 0.5]);
x = -b:0.01:a;     % array of the coordinate along beam length
y2 = scale_factor*(psi_2(1) - x*psi_2(2)); % vertical coordinates of the beam (2nd mode)

% plot mode shapes
plot(x,y2,'-k','linewidth',3);hold on  % 2nd mode
plot([x2_0], [0],'xr','linewidth',2,'markersize',12) % nodes

% Plot horizontal line
xlim1 = xlim; ylim1 = ylim;
line(xlim1.',[0; 0],...
    'linewidth',1,...
    'color',[0,0,0],'linestyle',':')
text([x(1), x(end), 0], [y2(1),y2(end),scale_factor*psi_1(1)],{'B','A','G'}, 'HorizontalAlignment','center',...
    'BackgroundColor',[1 1 1],'FontSize',12)
% Plot vertical line
axis equal
title(['2nd mode shape(f_2=',num2str(round(f_n(2)*100)/100),' Hz)']); xlabel('x [m]'); ylabel('y'); 

% Plot unsprung mass
hold on;
R = 0.1;
% offset Z_m_ns = -0.3;
h3 = rectangle('Position',[-b-R,offsetZ_m_ns-R,2*R,2*R],'FaceColor','none','Curvature',[0,0],'linestyle','--','linewidth',0.5);
h4 = rectangle('Position',[-b-R,offsetZ_m_ns+scale_factor*psi_2(4)-R,2*R,2*R],'FaceColor','b','Curvature',[0.8,0.4]);
h5 = rectangle('Position',[a-1*R,offsetZ_m_ns-R,2*R,2*R],'FaceColor','none','Curvature',[0,0],'linestyle','--','linewidth',0.5);
h6 = rectangle('Position',[a-1*R,offsetZ_m_ns+scale_factor*psi_2(3)-R,2*R,2*R],'FaceColor','b','Curvature',[0.8,0.4]);
grid on
% wheels
R = 0.3;
h3 = rectangle('Position',[-b-R,1*offsetZ_m_ns-R,2*R,2*R],'FaceColor','none','Curvature',[1,1]);
h5 = rectangle('Position',[a-1*R,1*offsetZ_m_ns-R,2*R,2*R],'FaceColor','none','Curvature',[1,1]);
% road
area(xlim,offsetZ_m_ns-R*[1 1]-0.3,offsetZ_m_ns-R,'linewidth',1,'FaceColor',[128 128 128]./256); hold on;

%% Mode 3
scale_factor = 2e-4;

h_fig_MS3 = figure('Name','Modal shapes','units','normalized','outerposition',[0 0 0.5 0.5]);
x = -b:0.01:a;     % array of the coordinate along beam length
y3 = scale_factor*(psi_3(1) - x*psi_3(2)); % vertical coordinates of the beam (1st mode)

% plot mode shapes
plot(x,y3,'-k','linewidth',3); hold on % 1st mode
plot([x3_0], [0],'xr','linewidth',2,'markersize',12) % nodes

% extension up to the node of mode 1
x_node = x3_0:0.01:-b;
y3_node = scale_factor*(psi_3(1) + x_node*psi_3(2));
plot(x_node,y3_node,'--b','linewidth',1); hold on

% Plot horizontal line
xlim1 = xlim; ylim1 = ylim;
line(xlim1.',[0; 0],...
    'linewidth',1,...
    'color',[0,0,0],'linestyle',':')
text([x(1),0, x(end)], [y3(1),scale_factor*psi_3(1),y3(end)],{'B','G','A'}, 'HorizontalAlignment','center',...
    'BackgroundColor',[1 1 1],'FontSize',12)
% Plot vertical line
axis equal
title(['3rd mode shape(f_3=',num2str(round(f_n(3)*100)/100),' Hz)']);xlabel('x [m]'); ylabel('y'); 

% Plot unsprung mass
hold on;
R = 0.1;
h3 = rectangle('Position',[-b-R,offsetZ_m_ns-R,2*R,2*R],'FaceColor','none','Curvature',[0,0],'linestyle','--','linewidth',0.5);
h4 = rectangle('Position',[-b-R,offsetZ_m_ns+scale_factor*psi_3(4)-R,2*R,2*R],'FaceColor','b','Curvature',[0.8,0.4]);
h5 = rectangle('Position',[a-1*R,offsetZ_m_ns-R,2*R,2*R],'FaceColor','none','Curvature',[0,0],'linestyle','--','linewidth',0.5);
h6 = rectangle('Position',[a-1*R,offsetZ_m_ns+scale_factor*psi_3(3)-R,2*R,2*R],'FaceColor','b','Curvature',[0.8,0.4]);
grid on
% wheels
R = 0.3;
h3 = rectangle('Position',[-b-R,1*offsetZ_m_ns-R,2*R,2*R],'FaceColor','none','Curvature',[1,1]);
h5 = rectangle('Position',[a-1*R,1*offsetZ_m_ns-R,2*R,2*R],'FaceColor','none','Curvature',[1,1]);
% road
area(xlim,offsetZ_m_ns-R*[1 1]-0.3,offsetZ_m_ns-R,'linewidth',1,'FaceColor',[128 128 128]./256); hold on;


%% Mode 4
scale_factor = 2e-4;

h_fig_MS4 = figure('Name','Modal shapes','units','normalized','outerposition',[0.5 0 0.5 0.5]);
x = -b:0.01:a;     % array of the coordinate along beam length
y4 = scale_factor*(psi_4(1) - x*psi_4(2)); % vertical coordinates of the beam (2nd mode)

% plot mode shapes
plot(x,y4,'-k','linewidth',3);hold on  % 2nd mode
plot([x4_0], [0],'xr','linewidth',2,'markersize',12) % nodes

% Plot horizontal line
xlim1 = xlim; ylim1 = ylim;
line(xlim1.',[0; 0],...
    'linewidth',1,...
    'color',[0,0,0],'linestyle',':')
text([x(1), x(end), 0], [y4(1),y4(end),scale_factor*psi_4(1)],{'B','A','G'}, 'HorizontalAlignment','center',...
    'BackgroundColor',[1 1 1],'FontSize',12)
% Plot vertical line
axis equal
title(['4th mode shape(f_4=',num2str(round(f_n(4)*100)/100),' Hz)']); xlabel('x [m]'); ylabel('y'); 

% Plot unsprung mass
hold on;
R = 0.1;
h3 = rectangle('Position',[-b-R,offsetZ_m_ns-R,2*R,2*R],'FaceColor','none','Curvature',[0,0],'linestyle','--','linewidth',0.5);
h4 = rectangle('Position',[-b-R,offsetZ_m_ns+scale_factor*psi_4(4)-R,2*R,2*R],'FaceColor','b','Curvature',[0.8,0.4]);
h5 = rectangle('Position',[a-1*R,offsetZ_m_ns-R,2*R,2*R],'FaceColor','none','Curvature',[0,0],'linestyle','--','linewidth',0.5);
h6 = rectangle('Position',[a-1*R,offsetZ_m_ns+scale_factor*psi_4(3)-R,2*R,2*R],'FaceColor','b','Curvature',[0.8,0.4]);
grid on
% wheels
R = 0.3;
h3 = rectangle('Position',[-b-R,1*offsetZ_m_ns-R,2*R,2*R],'FaceColor','none','Curvature',[1,1]);
h5 = rectangle('Position',[a-1*R,1*offsetZ_m_ns-R,2*R,2*R],'FaceColor','none','Curvature',[1,1]);
% road
area(xlim,offsetZ_m_ns-R*[1 1]-0.3,offsetZ_m_ns-R,'linewidth',1,'FaceColor',[128 128 128]./256); hold on;

%% Modal analysis of the damped system
% eigenvalues and eigenvectors (both complex)
n = 4;   % number of system degrees of freedom
A = [ zeros(n,n)  diag(ones(n,1))
    -mm^-1*kk       -mm^-1*cc]; % state matrix of state-space representation
[V,D] = eig(A)  % produces a diagonal matrix lambda of generalized eigenvalues and a
%     full matrix vv whose columns are the corresponding eigenvectors so that A*V = V*D
[Wn,Z,P] = damp(A); 
% [Wn,Z,P] = damp(SYS) returns vectors Wn, Z and P containing the natural frequencies, damping factors and poles of the linear system SYS.
f_damp = Wn/(2*pi)
abs(P)/(2*pi)
Z

%% Sensitivity analysis of the damping ratio for different viscous coefficient
flag_plot_comp_mode = 0;

switch flag_plot_comp_mode
    case 1
        h_fig_mode_comp = figure('units','normalized','Outerposition',[0 0 0.5 1]); hold all
        h_fig_freq = figure('units','normalized','Outerposition',[0.5 0 0.5 1]); hold all
end

cs1_vet = cs1*linspace(0.01,1.5,10000);
cs2_vet = cs2*linspace(0.01,1.5,10000);

cs1_n = cs1;
cs2_n = cs2;

% Initialization
f_c = zeros(length(cs1_vet),length(P));
zeta = zeros(length(cs1_vet),length(P));

for cont1=1:length(cs1_vet)
    cs1 = cs1_vet(cont1); cs2 = cs2_vet(cont1);
    
    % Damping matrix update
    cc = [cs1+cs2   -cs1*a+cs2*b         -cs1    -cs2
        -cs1*a+cs2*b   cs1*a^2+cs2*b^2   cs1*a   -cs2*b
        -cs1               cs1*a          c1+cs1     0
        -cs2               -cs2*b             0       c2+cs2];

    A = [ zeros(n,n)  diag(ones(n,1))
        -mm^-1*kk       -mm^-1*cc];
    [V,D] = eig(A);  % produces a diagonal matrix lambda of generalized eigenvalues and a
    %     full matrix vv whose columns are the corresponding eigenvectors so that A*V = V*D
   
    psi_1 = V(:,1)./(V(1,1));   % 1st eigenvector
    psi_2 = V(:,2)./(V(1,2));   % 2nd eigenvector
    psi_3 = V(:,3)./(V(1,3));   % 1st eigenvector
    psi_4 = V(:,4)./(V(1,4));   % 2nd eigenvector
    
    psi = [psi_1 psi_2 psi_3 psi_4];
    
    [Wn,Z,P] = damp(A);
    f_damp = Wn/(2*pi);
    
    f_c(cont1,:) = [f_damp];
    zeta(cont1,:) = [Z];
end

flag_sort_f = 0;
switch flag_sort_f
    case 1
        [f_c_ord, in] = sort(f_c,2);
        for cont2=1:length(in)
            in(cont2,:);
            zeta_ord(cont2,:) = zeta(cont2,in(cont2,:));
        end
    case 0
        f_c_ord = f_c;
        zeta_ord = zeta;
end

%% Plot f e zeta(c)
figure('units','normalized','outerposition',[0 0.5 0.5 0.5]); 
plot(f_c_ord(:,2).*sqrt(1-zeta_ord(:,2).^2),cs1_vet/cs1_n,'or','linewidth',2); hold on
plot(f_c_ord(:,4).*sqrt(1-zeta_ord(:,4).^2),cs1_vet/cs1_n,'.k','linewidth',3);
plot(f_c_ord(:,6).*sqrt(1-zeta_ord(:,6).^2),cs1_vet/cs1_n,'+b','linewidth',4);
plot(f_c_ord(:,8).*sqrt(1-zeta_ord(:,8).^2),cs1_vet/cs1_n,'dm','linewidth',1);
plot(xlim,[1 1],'--k')
% ylim([0 5])
leg = legend('f_2','f_4','f_6','f_8','Location','best'); set(leg,'Fontsize',14)
ylabel('normalized dampers viscous friction coefficient: c_s/c_{s,nom}','Fontsize',14); xlabel('Damped natural frequency: f_d [Hz]'); grid on;
set(gca,'Fontsize',14)

% Plot zeta
figure('units','normalized','outerposition',[0.5 0.5 0.5 0.5]); 
plot(zeta_ord(:,2)*100,cs1_vet/cs1_n,'or','linewidth',1); hold on
plot(zeta_ord(:,4)*100,cs1_vet/cs1_n,'.k','linewidth',1); 
plot(zeta_ord(:,6)*100,cs1_vet/cs1_n,'+b','linewidth',1); 
plot(zeta_ord(:,8)*100,cs1_vet/cs1_n,'dm','linewidth',1); 
plot(xlim,[1 1],'--k')
ylabel('normalized dampers viscous friction coefficient: c_s/c_{s,nom}','Fontsize',14); xlabel('Damping factor: \zeta [%]'); grid on; 
leg = legend('\zeta_2','\zeta_4','\zeta_6','\zeta_8','Location','best'); set(leg,'Fontsize',14)
set(gca,'Fontsize',14)

%% FRF 
% Damping
cc = [cs1_n+cs2_n   -cs1_n*a+cs2_n*b         -cs1_n    -cs2_n
    -cs1_n*a+cs2_n*b   cs1_n*a^2+cs2_n*b^2   cs1_n*a    -cs2_n*b
    -cs1_n                cs1_n*a          c1+cs1_n     0
    -cs2_n               -cs2_n*b             0       c2+cs2_n];

% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% System response (INVERSE of the dynamic stiffness matrix)  
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%Variables and initialization
V = 100/3.6;
tau = L/V;

OM = 2*pi*[0:1e-2:30]';      % rad/s
lambda = 2*pi*V./OM;         % m

Ku = [0 0; 0 0; k1 0; 0 k2];
Cu = [0 0; 0 0; c1 0; 0 c2];

%  initialization of FRF vectors
Z_G_ddot      = ones(length(OM),1);
theta_ddot    = ones(length(OM),1);
Z_F_ddot      = ones(length(OM),1);
Z_R_ddot      = ones(length(OM),1);

%%%%%%%%%%%%%%%%%%% MY PROCEDURE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h0 = 1/100; %1 cm of amplitude, seen in the slides

%% SYNCRONOUS
tau = 0;
for count1=1:length(OM)
    ww      = OM(count1);

    u1     = 1;
    u2     = 1;    %To obtian the trasnfer function of the system of our 
    % interest it is necessary to consider the input that has been imposed.
    % In syncronous case it is u0 defined as unit column
    Kdyn    = kk + i*ww*cc - ww^2*mm;
    Rec     = inv(Kdyn);
    Trans   = Rec*(Ku + i*ww*Cu);
    
    Z_G_ddot(count1,1) = -ww^2*Trans(1,:)*[u1 ; u2];
    theta_ddot(count1,1) = -ww^2*Trans(2,:)*[u1 ; u2];
    Z_F_ddot(count1,1)= -ww^2*Trans(3,:)*[u1 ; u2];
    Z_R_ddot(count1,1)= -ww^2*Trans(4,:)*[u1 ; u2];
end

% Plots

figure(7); 
subplot(221); hold all;
plot(OM/(2*pi),abs(Z_G_ddot),'Linewidth',1.5); grid on;

subplot(222); hold all;
plot(OM/(2*pi),abs(theta_ddot),'Linewidth',1.5); grid on;

subplot(223); hold all;
plot(OM/(2*pi),abs(Z_F_ddot),'Linewidth',1.5); grid on;

subplot(224); hold all;
plot(OM/(2*pi),abs(Z_R_ddot),'Linewidth',1.5); grid on;

c = 4.7*10^-4; n=2.1;
f = OM./(2*pi);

S = c*V^(n-1).*f.^-n;
S_r = Z_G_ddot.^2.*S;

figure(9); loglog(f,abs(S)); grid on; hold on;
figure(10); loglog(f,abs(S_r)); grid on; hold on;
figure(11); plot(f,abs(S_r)); grid on; hold on;

%% DELAYED

Z_G_ddot      = ones(length(OM),1);
theta_ddot    = ones(length(OM),1);
Z_F_ddot      = ones(length(OM),1);
Z_R_ddot      = ones(length(OM),1);
tau = L/V;

for count1=1:length(OM)
    ww      = OM(count1);

    u1     = 1;
    u2     = exp(-i*ww*tau);
    
    Kdyn    = kk + i*ww*cc - ww^2*mm;
    Rec     = inv(Kdyn);
    Trans   = Rec*(Ku + i*ww*Cu);
    
    Z_G_ddot(count1,1) = -ww^2*Trans(1,:)*[u1 ; u2];
    theta_ddot(count1,1) = -ww^2*Trans(2,:)*[u1 ; u2];
    Z_F_ddot(count1,1)= -ww^2*Trans(3,:)*[u1 ; u2];
    Z_R_ddot(count1,1)= -ww^2*Trans(4,:)*[u1 ; u2];
end

% Plots
figure(7); 
subplot(221); hold all;
plot(OM/(2*pi),abs(Z_G_ddot),'r--','Linewidth',1.5); grid on;
ylabel('a_{zg}\h [(m/s^2)/m]')

subplot(222); hold all;
plot(OM/(2*pi),abs(theta_ddot),'r--','Linewidth',1.5); grid on;
ylabel('\theta_{dd}/h[(rad/s^2)7m]')

subplot(223); hold all;
plot(OM/(2*pi),abs(Z_F_ddot),'r--','Linewidth',1.5); grid on;
xlabel('Frequency [Hz]')
ylabel('a_{z1}/h [(m/s^2)/m]')

subplot(224); hold all;
plot(OM/(2*pi),abs(Z_R_ddot),'r--','Linewidth',1.5); grid on;
xlabel('Frequency [Hz]')
ylabel('a_{z2}/h [(m/s^2)/m]')
legend('Syncronous','Delayed','Location','best')

c = 4.7*10^-4; n=2.1;
f = OM./(2*pi);

S = c*V^(n-1).*f.^-n;
S_r = Z_G_ddot.^2.*S;

figure(9); loglog(f,abs(S)); grid on; legend('Syncronous','Delayed','Location','best')
figure(10); loglog(f,abs(S_r)); grid on; legend('Syncronous','Delayed','Location','best')
figure(11); plot(f,abs(S_r)); grid on; legend('Syncronous','Delayed','Location','best')

%% PSD without for
% c = 4.7*10^-4; n=2.1;
% f = OM./(2*pi);
% 
% S = c*V^(n-1).*f.^-n;
% S_r = Z_G_ddot.^2.*S;
% 
% figure(9); loglog(f,abs(S)); grid on; legend('S')
% figure(10); loglog(f,abs(S_r)); grid on; legend('S_r')
% figure(11); plot(f,abs(S_r)); grid on; legend('S_r no loglog')

%% PSD with for cycle
% S = ones(length(OM),1);
% S_r = ones(length(OM),1);
% f = OM./(2*pi);
% c = 4.7*10^-4; n=2.1;
% 
% for count1=1:length(OM)
%     w = OM(count1,1);
%     S(count1) = c*V^(n-1)*(w/(pi*2))^-n;
%     S_r(count1) = Z_G_ddot(count1)^2*S(count1);
% end
% 
% figure(9); loglog(f,abs(S)); grid on; legend('S')
% figure(10); loglog(f,abs(S_r)); grid on; legend('S_r')
% figure(11); plot(f,abs(S_r)); grid on; legend('S_r no loglog')