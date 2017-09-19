function CVcongest(dx_trans, plotArgs);
% CVcongest - congestion of traveling wave

[dx_trans, plotArgs] = arginDefaults('dx_trans,plotArgs',1,'n');

NN = numel(dx_trans);
if NN>1, % recursion
    plotArgs = ['n' ploco(2:NN)];
    for ii=1:NN,
        dx = dx_trans(ii);
        CVcongest(dx, plotArgs(ii));
    end
    subplot(2,2,2);
    legend(Words2cell(sprintf( '%0.1f$',  dx_trans),'$'), 'location', 'northwest');
    return;
end

% =========single dx_trans from here==========

U1 = 10.01; % m/s fast wave 
U2 = 1; % m/s slow wave
xA = 4;
xB = 6;
Fmin = 5;
Fmax = 15;

if isequal('n', plotArgs), % new fig
    set(figure,'units', 'normalized', 'position', [0.184 0.265 0.46 0.642]);
    plotArgs='b';
end

x = linspace(0, 10).';
x_trans = linspace(1,9,150);
Freq = interp1(x_trans([end 1]), [Fmin Fmax], x_trans);
U = local_U(x,x_trans, dx_trans, U1, U2);
[dum, iA] = min(abs(x-xA)); % index of xA in x
[dum, iB] = min(abs(x-xB)); % ditto xB
xmid = (xA+xB)/2;
[dum, ix_mid] = min(abs(x_trans-xmid));

subplot(2,2,1);
xplot(x, U(:,[1 ix_mid end]), plotArgs);
ylim([0, 1.1*U1]);
fenceplot([xA xB], ylim, 'color', 0.8*[1 1 1], 'linewidth', 3);
xlabel('BM location');
ylabel('Group velocity');

[tau, dx] = local_tau(xA,xB,x,U);
subplot(2,2,3);
xplot(Freq,tau, plotArgs);
xlabel('Frequency');
ylabel('Group delay A->B');

subplot(2,2,2);
Gain = P2dB(U(iA,:)./U(iB,:));
xplot(Freq,Gain, plotArgs);
xlabel('Frequency');
ylabel('Gain (dB)');

subplot(2,2,4);
df = abs(mean(diff(Freq)));
PH = df*cumsum(tau);
PH = PH-PH(end);
xplot(Freq, PH, plotArgs);
phasedelayslider(gca,'kHz');
xlabel('Frequency');
ylabel('Phase A->B');


function U = local_U(x, x_trans, dx_trans, U1, U2);
[x, x_trans] = SameSize(x, x_trans);
x1 = x_trans - dx_trans/2; % start of transition
x2 = x_trans + dx_trans/2; % end of transition
%[U1 U2] = dealelements(log([U1 U2]));
U = U1*double(x<=x1) + U2*double(x>=x2) ...
    + (U1+((x-x1).*(U2-U1)./(x2-x1))).*betwixt(x,x1,x2);
%U = (U1+U2 + (U2-U1)*2*atan(5*(x-x_trans)./dx_trans)/pi)/2;
%U = exp(U);

function [tau, dx] = local_tau(xA,xB,x,U);
% travel time from A to B
dx = diff(x(1:2));
U(x<xA,:) = inf; % this yields zero travel time before A
U(x>xB,:) = inf; % this yields zero travel time beyond B
%plot(x,dx./U);
tau = dx*sum(1./U,1);

