function [freq, A] = artOHC(Cm, Rm, Ra, plotArg);
% artOHC - complex admittance of model OHC 
%    artOHC(Cm, Rm, Ra) plots complex admittance (I-vs-V transfer) of a
%    linear OHC model OHC consisting of an access resistance Ra in series
%    with a membrane having resistance Rm and capacitance Cm
%    Rm, Ra in Mohm, Cm in pF.
%    Defaults:
%       Cm = 60 pF
%       Rm = 100 Mohm
%       Ra = 5 Mohm
%
%    See also ...

if nargin<1, Cm = 60; end;
if nargin<2, Rm = 100; end;
if nargin<3, Ra = 5; end;
if nargin<4, plotArg='n'; end;
[Cm, Rm, Ra] = SameSize(Cm, Rm, Ra);


%=======single param set from here============
freq = logispace(5,10e3,1e4);
Am = 1e-6/Rm + 1e-12*Cm*i*2*pi*freq; % admittance of membrane 
A = impara(Am,1e-6/Ra); % Ra in series
I = 1e9*0.1*A; % current [nA] in response to a 100-mV stimulus.

subplot(2,1,1);
xplot(freq,A2dB(abs(I)), plotArg);
xlog125;
subplot(2,1,2);
xplot(freq,cunwrap(angle(I)/2/pi), plotArg);
xlog125;







