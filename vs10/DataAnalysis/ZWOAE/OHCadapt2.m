function [T, MET, S] = OHCadapt2(tau, A, freq, plotArg);
% OHCadapt - toy model for OHC MET channel adaptation
%    [Stim Resp] = OHCadapt(tau, A, freq);
%    Inputs
%        tau: [us] adaptation time constant
%          A: [AU] stimulus amplitude
%       freq: [Hz] stimulus freq
%    Outputs
%       Stim: one cycle of the cosine stimulus as specified
%       Resp: one cycle of the steady state response

%IOfcn = @(x)max(-2, min(x,1));
IOfcn = @tanh;
%IOfcn = @(x)x;;
%IOfcn = @(x)(0.5+tanh(x)).^2;
% IOfcn = @(x)(0.5*(1+tanh(x+0.5))).^3;

if nargin<4, plotArg=''; end
if isempty(plotArg), initplot = @cla; else, initplot=@nope; end

% initialize local_dydt
S.tau = tau;
S.Stim = @(t)A*sin(2e-3*pi*t*freq);
%S.Stim = @(t)A*double((t>10) & (t<20));
local_dydt(S); % initialization
% solve ODE
Ncyc = 20;
T = linspace(0, Ncyc*1e3/freq, Ncyc*50); % 
A = 0*T;
Stim = IOfcn(S.Stim(T));
dt = diff(T(1:2))
for ii=1:numel(T)-1,
    dy = dt*local_dydt(A(ii), Stim(ii));
    A(ii+1) = A(ii)+ dy;
end

%MET = IOfcn(Stim-A/2);
MET = (Stim-A/2);
f1; set(gcf,'units', 'normalized', 'position', [0.0906 0.511 0.438 0.41]);
plot(T,Stim,'-r')
xplot(T,MET);
xplot(T,A,'-g')
f2; set(gcf,'units', 'normalized', 'position', [0.541 0.51 0.438 0.41])
plot(Stim, MET);
%xplot(T,A,':g')
%xplot(T,local_dydt(T,A),'k')
%=======================
function dy = local_dydt(y,x);
persistent S
if isstruct(y), % initialize
    S = y;
    S.mtau = 1e-3*mean(S.tau);
    S.dtau = -1e-3*0.5*diff(S.tau);
    S
    return
end
D = (x-y);
tau = S.mtau + sign(D)*S.dtau;
dy = D./tau;
%xplot(t,x,'.')


