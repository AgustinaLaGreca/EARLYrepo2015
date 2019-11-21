function report_alt_level(figh, Prefix, NS)
% report_alt_level - report spectrum or total SPL, whichever is not specified
%   report_alt_level(figh, Prefix, SPLUnit, NS)
%  Inputs
%       figh: handle to stim menu
%      Prefix: name prefix of the noise params (e.g. Mnoise fro MASK stim)
%        NS: output of NoiseSpec as used by the stimulus generator (makestimXXX)
%  report_alt_level reports to the 

if ~ishandle(figh), return; end
M = find(messenger(), figh, [Prefix, 'AltLevel']);
if isempty(M), return; end
totSPL = unique([NS.totSPL]);
tStr = trimspace(num2str(0.1*round(totSPL*10)))
specSPL = unique([NS.specSPL]);
sStr = trimspace(num2str(0.1*round(specSPL*10)))
switch NS(1).SPLUnit,
    case 'dB/Hz',
        %Str = sprintf('', );
    case 'dB SPL',
end


