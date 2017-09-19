function S = coch3Dk(RR, CC, B, plotArgs);
% coch3Dk - frequency versus wavenumber, 3D wave analysis
%   coch3Dk(RR, CC, plotArgs);

[RR, CC, B, plotArgs] = arginDefaults('RR, CC, B, plotArgs',0, 1, [], []);

L = 0.5; % mm scale radius
if isempty(B),
    B = logispace(5, 100, 5); % um effective BM width
end
R = 1e-3*B/pi; % mm effective BM "radius" (width / pi)
lambda = logispace(0.01, 100).'; % mm wave length of traveling wave
k = 2*pi./lambda; % wavenumber 
om = logispace(0.1,1e3,1e4).'; % angular freq 

[RR, CC, B, R] = SameSize(RR, CC, B, R);

if isempty(plotArgs),
    figure;
    plotArgs = struct('linestyle', '-');
    LegStr = {};
end
set(gcf,'units', 'normalized', 'position', [0.0805 0.175 0.863 0.71]);
LegStr = getGUIdata(gcf, 'LegStr', {});
for ii=1:numel(B),
    r = R(ii); % radius
    cc = CC(ii); % compliance
    rr = RR(ii); % resistance
    LegStr{end+1} = sprintf('RR=%0.2f  C = %0.1f  B=%d', rr, cc, round(B(ii)));
    M(:,ii) = 0.01*k.*SteeleForm(k*L,k*r)/r;
    OM1(:,ii) = k./sqrt(M(:,ii)*cc);
    %OM2(:,ii) = k.^2./M(:,ii);
    Gom = om.*real((1+i*om*rr).^(-0.5));
    Gang = cunwrap(cangle((1+i*om*rr).^(-0.5)));
    OM2(:,ii) = interp1(Gom, om, OM1(:,ii));
    PH_k(:,ii) = interp1(Gom, Gang, OM1(:,ii));
    U2(:,ii) = diff(OM2(:,ii))./diff(k);
end
setGUIdata(gcf, 'LegStr', LegStr);
%=
k = SameSize(k, M);
subplot(2,3,1);
xplot(k, OM2, plotArgs);
xlog125([0.1 100]);
ylog125([0.1 100]);
xlabel('k');
ylabel('\omega');
legend(LegStr, 'location', 'northwest')
%=
subplot(2,3,2);
xplot(OM2(2:end,:), P2dB(U2), plotArgs);
xlog125([0.1 100]);
grid on
xlabel('\omega');
ylabel('U (dB)');
%=
subplot(2,3,3);
xplot(OM2(2:end,:), diff(log(OM2))./diff(log(k)), plotArgs)
xlog125([0.1 100]);
ylog125([0.1 5]);
grid on
xlabel('\omega');
ylabel('U/c');
%=
subplot(2,3,4);
xplot(k(2:end, :), P2dB(U2), plotArgs)
xlog125([0.1 100]);
grid on
xlabel('k');
ylabel('U (dB)');
%=
subplot(2,3,5);
xplot(OM2(2:end-1, :), diff(log(U2))./diff(log(OM2(2:end,:))), plotArgs)
xlog125([0.1 100]);
grid on
xlabel('\omega');
ylabel('dlog(U)/dlog(\omega)');
%=
subplot(2,3,6);
xplot(OM2, PH_k, plotArgs)
xlog125([0.1 100]);
set(gca, 'ytick', 1/32*(-5:1), 'ylim', [-5 1]/32);
grid on
xlabel('\omega');
ylabel('angle(k)');
%


S = CollectInStruct('-primary', L, B, R, lambda, k, '-derived', M, OM1, OM2, U2);

