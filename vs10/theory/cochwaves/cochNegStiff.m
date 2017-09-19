% ===sandbox for 3D neg st model=====

k = logispace(1e-3,1,1e5); 
SF = steeleForm(k*10,k/2); % 3D mass factor a la Steele & Taber
% reference: no stiffness
U0=0; alpha = 0; beta = 0; k(k>0.51)=nan; omega=sqrt(k./SF.*(1-alpha*(sin(pi*k).^4-beta*(pi*k).^4))); U = diff(omega)./diff(k); c = omega./k; U0=U; omega0 = omega;

f1; plot(omega0(2:end),U0,'k'); xlog125; grid on; xlim(max(omega0)*[0.1 2]); ylim([-5 20])

% markers: imin(3) is wavenumber  corresponding to 8*longitudinal coupling
% via Deiters' phalanxes
[dum imin(1)] = min(abs(k-0.5)); [dum imin(2)] = min(abs(k-0.25)); [dum imin(3)] = min(abs(k-0.125));

alpha = 50; beta = 0.8; omega=sqrt(k./SF.*(1-alpha*(sin(pi*k).^4-beta*(pi*k).^4))); U = diff(omega)./diff(k); c = omega./k;
xplot(omega(2:end),U,'color', rand(1,3)); xplot(omega(imin+1), U(imin),'r*'); xplot(omega(imin(end)+1), U(imin(end)),'ks');
alpha = 100; beta = 0.82; omega=sqrt(k./SF.*(1-alpha*(sin(pi*k).^4-beta*(pi*k).^4))); U = diff(omega)./diff(k); c = omega./k;
xplot(omega(2:end),U,'color', rand(1,3)); xplot(omega(imin+1), U(imin),'r*'); xplot(omega(imin(end)+1), U(imin(end)),'ks');

% various plots
f2; plot(omega(2:end),U./c(2:end)); xlog125; ylim([0 3]); min(U./c(2:end))
f3; plot(omega,c)
f4; plot(k,omega)
f5; plot(omega,-k); xlim([0 1])
f6; plot(k(2:end), U); xlim([0 0.3]);







