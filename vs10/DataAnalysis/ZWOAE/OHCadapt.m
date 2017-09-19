function [Stim, Resp] = OHCadapt(tau, A, freq, plotArg);
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
%IOfcn = @(x)(0.5+tanh(x)).^2;
% IOfcn = @(x)(0.5*(1+tanh(x+0.5))).^3;

if nargin<4, plotArg=''; end
if isempty(plotArg), initplot = @cla; else, initplot=@nope; end

[tau, A, freq] = SameSize(tau, A, freq);
if numel(tau)>1, % recursion
    [Stim, Resp] = deal([]);
    for ii=1:numel(tau),
        [stim, resp] = OHCadapt(tau(ii), A(ii), freq(ii), ploco(ii));
    end
    Stim = [Stim, stim];
    Resp = [Resp, resp];
    return;
end

Nsam = 1000;
phi = 2*pi*(0:Nsam-1).'/Nsam; % phase vector [rad]
Zstim = -i*A*exp(i*phi); % one stimulus cycle
i0 = 1;%round(Nsam/4);

om_tau = 2*pi*freq*1e-6*tau; % omega*tau in cycles
DX = Zstim*(1+i*om_tau/(1+i*om_tau))/2; % complex Stimulus re adapting resting position
Resp = IOfcn(real(DX));
Stim = real(Zstim);
% f3; set(gcf,'units', 'normalized', 'position', [0.436 0.212 0.44 0.28])
% xplot(phi, Stim, plotArg);

f2; set(gcf,'units', 'normalized', 'position', [0.438 0.595 0.438 0.31]);
xplot(phi, Resp, plotArg);

f1;set(gcf,'units', 'normalized', 'position', [0.0148 0.351 0.367 0.551]);
initplot();
xplot(Stim,Resp, plotArg);
xplot(Stim(i0),Resp(i0),plotArg, 'marker', 'o');
xplot(Stim(i0+15),Resp(i0+15),plotArg, 'marker', '*');

% f4;
ph1 = angle(Zstim'*Resp); % hysteresis angle 
% xplot(real(Zstim*exp(i*ph1)), Resp, plotArg);
% f4;
% xplot(phi/2/pi,Stim, 'b');
% xplot(phi/2/pi,Stim-real(DX),'r');
% xplot(phi/2/pi,real(DX),'g');




