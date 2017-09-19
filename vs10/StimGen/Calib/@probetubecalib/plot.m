function ah=plot(P, varargin);
% Probetubecalib/plot - plot Probetubecalib object
%   plot(PTC) plots ProbeTubeCalib object PTC, displaying the Magnitude and
%   phase transfers of the Probe & Cavity recordings in one pair of graphs
%   and those of the  Probe-to-Cavity transfer in aonther pair. Plot
%   returns graphics handles to the four axes systems.
%
%   See also Transfer/Plot.

% equalize the WB dleays of Probe & Cavity
setWBdelay(P.Probe, getWBdelay(P.Cavity));
% create subpolots
ah = [subplot(2,2,1),subplot(2,2,3),subplot(2,2,2),subplot(2,2,4)];
figh = gcf;
set(figh,'units', 'normalized', 'position', [0.324 0.286 0.666 0.628]);
% plot Probe & Cavity together after synching their delay
WBdelay = max(getWBdelay([P.Probe P.Cavity]));
P.Probe = setWBdelay(P.Probe,WBdelay);
P.Cavity = setWBdelay(P.Cavity,WBdelay);
hP=plot(ah(1:2), P.Probe, varargin{:});
hC=plot(ah(1:2), P.Cavity, 'r', varargin{:});
legend([hP(1) hC(2)], 'Probe', 'Cavity');
% plot Tube
plot(ah(3:4),P.Tube, varargin{:});
title(ah(3), P.Descr, 'interpreter', 'none');




