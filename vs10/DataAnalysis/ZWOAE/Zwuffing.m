function Mess = ZWOAEtest(Ls, zw,PlotArg);
% ZWuffing - Suppression by a Duffing oscillator
%   ZWuffing(Ls, [Fzw [PlotArg]], 'plotarg')

if nargin<1, Ls=40; end  % default 40 dB single tone
if nargin<2, zw=[0 -inf]; end  % default: no zwuis group
if nargin<3, PlotArg='n'; end

tspan = [0 2010];
y0 = [0; 0];

Fzw = zw(1);
Lzw = zw(2);

DF = cumsum(10+[0:6])/1000;
Nzw = numel(DF);
DF = DF-DF(round(Nzw/2));
Fzw= 0.001*round(1000*(Fzw+DF)); % rel freq of zwuis group re single tone
om = 2*pi; % angular resonance freq
omx = om*[1 Fzw]; % stimulus: single tone (freq=1) and zwuis below it
L = [Ls Lzw*ones(1,Nzw)].'; 
Ampl = dB2A(L);
nn = 0.1; % nonlinearity index
alpha = 0.1; % friction coeff
dydt = @(t,y)[y(2); cos(omx*t)*Ampl-alpha*om*y(2)-om^2*(y(1)+nn*y(1).^3)];
% solve the problem using ODE45
N = 2e4;
[t,y]=ode45(dydt,tspan,y0);
%dsize(t,y)
%plot(t,y);
Dtrans = 10; % transient portion of response to be tossed
Tsteady = linspace(Dtrans, tspan(2), N+1).';
Tsteady = Tsteady(1:end-1); % commensurate w stim freqs
Y = interp1(t,y(:,1),Tsteady);
SNR=10;
figure(1);
set(gcf,'units', 'normalized', 'position', [0.091875 0.49583 0.43938 0.41167])
X = cos(Tsteady(:)*omx)*Ampl;
plot(Tsteady,X/std(X)); xplot(Tsteady,Y/std(Y),'r');legend('stim','response'); 
Y = Y+randn(size(Y))*dB2A(-SNR)*std(Y);
figure(2);
set(gcf,'units', 'normalized', 'position', [0.54563 0.49 0.43938 0.41167])

subplot(2,1,1);
% plot(Tsteady,Y);
% magn spectrum
df = 1/(max(Tsteady)-min(Tsteady));
FreqAxis = df*(0:N-1);
Sp = fft(Y(:,1));
MG = A2dB(abs(Sp)); 
PH = angle(Sp)/2/pi; 
PH = delayPhase(PH, FreqAxis(:), -Dtrans); % compensate timing shift due to deletion of transient portion
xdplot(df,MG,PlotArg);
ylim([-50 0]+max(MG));
xlim([0,4*om/2/pi]);

Fprim = [Fzw(:); 1]; % stimulus freqs
GDprim = local_GD(df,MG,PH, Fprim(1:end-1), ZWOAEplotStruct('Prim'))
local_GD(df,MG,PH, Fprim(end), ZWOAEplotStruct('Prim'))

Msup=ZWOAEmatrices(Nzw,'Msupall');
Fsup = Msup*Fprim; % freqs of supall DPs
GDsup = local_GD(df,MG,PH, Fsup, ZWOAEplotStruct('DPsupall'))

Ffar = ZWOAEmatrices(Nzw,'Mfar')*Fprim; % freqs of far DPs
GDfar = local_GD(df,MG,PH, Ffar, ZWOAEplotStruct('DPfar'))
% --------------------------------------------------------------------------

function GD=local_GD(df,MG,PH, Freq,PlotArg);
N = numel(MG); FreqAxis = df*(0:N-1);
isup = interp1(FreqAxis,1:numel(FreqAxis), Freq(:), 'nearest'); % indices of Freq comps
subplot(2,1,1);
xplot(Freq,MG(isup), PlotArg);
subplot(2,1,2);
ph = ucunwrap(PH(isup),Freq, '-tozero');
xplot(Freq, ph, PlotArg, 'linestyle', '-');
GD = polyfit(Freq,ph,1); GD = -GD(1);
text(min(Freq),0.1*(diff(ylim))+max(ph),['GD=' num2str(GD)], 'verticalalign','bottom');



