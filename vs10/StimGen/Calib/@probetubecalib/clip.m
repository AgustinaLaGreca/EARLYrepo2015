function P = Clip(P, Freq, MinVal, MaxVal);
% Probetubecalib/Clip - clip part of Tube response between specified limits
%   PTC = Clip(PTC, [f1 f2], Mn Mx) clips the magnitude of a restricted 
%   frequency range of the Tube trf of PTC. The values are limited 
%   between minimum value Mn dB and max value mx dB.
%   The frequency range is indicated by the frequencies f1 and f2 (Hz).
%   The change induced by clip is added to the EditHistory of PTC, unless
%   the action was void, e.g. an empty frquency range.
%
%   See also Probetubecalib/edit, Probetubecalib/smooth.

if ischar(MinVal), MinVal = eval(MinVal); end
if ischar(MaxVal), MaxVal = eval(MaxVal); end

ihit = find(betwixt(frequency(P.Tube), Freq(1), Freq(2)));
if isempty(ihit), return; end

MG = magnitude(P.Tube);
PH = phase(P.Tube);
MG(ihit) = clip(MG(ihit), MinVal, MaxVal);

% now put this back in P.Tube
tub = struct(P.Tube); % ugly, but fro now faster than writing a transfer/clip
tub.Ztrf = dB2A(MG).*exp(2*pi*i*PH);
P.Tube = transfer(tub);

P.EditHistory = [P.EditHistory struct('fcn', fhandle(mfilename), 'args', {{probetubecalib, Freq, MinVal, MaxVal}})];













