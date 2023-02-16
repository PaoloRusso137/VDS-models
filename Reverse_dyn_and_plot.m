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
plot(OM/(2*pi),abs(X12),'-.','DisplayName',sprintf('%d',c));
ylabel('X_{12} = z_s/h)'); legend('show'); grid on;

subplot(212); hold all; 
plot(OM/(2*pi),abs(X22),'-','DisplayName',sprintf('%d',c));
ylabel('X_{22} = z_u/h'); xlabel('frequency [Hz]')
legend('show'); grid on;

figure(2); 
subplot(211); hold all;
plot(OM/(2*pi),OM.^2.*abs(X12),':','DisplayName',sprintf('%d',c));
ylabel('zpp_s/h'); legend('show'); grid on;

subplot(212); hold all;
plot(OM/(2*pi),abs(Fz)./k1p,'--','DisplayName',sprintf('%d',c));
ylabel('F_z/(k_p h)'); xlabel('frequency [Hz]')
legend('show'); grid on;