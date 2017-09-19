function [Gain, Phase] = ProbeLoss(FN, Freq);
% ProbeLoss - transfer function of probe tube
%   [Gain, Phase] = ProbeLoss(ProbeName, Freq)
%   returns Gain [dB] and Phase [cycle] of probe-tube transfer function of
%   probe ProbeName, whose default extension is 'probeloss'. 
%   Freq is array of frequency [Hz].
%
%   See also ProbeTubeCalib, MeasureTransfer.

persistent last_FN last_TRF
if ~isequal(FN, last_FN),
    last_TRF = local_read(FN);
    last_FN = FN;
end

% find(isnan(last_TRF.Freq))
% find(isnan(last_TRF.Magn))
% find(isnan(last_TRF.Phase))
Gain = interp1(last_TRF.Freq, last_TRF.Magn, Freq);
Phase = interp1(last_TRF.Freq, last_TRF.Phase, Freq);
Gain(Freq<last_TRF.Freq(1)) = last_TRF.Magn(1);
Phase(Freq<last_TRF.Freq(1)) = last_TRF.Phase(1);
Gain(Freq>last_TRF.Freq(end)) = last_TRF.Magn(end);
Phase(Freq>last_TRF.Freq(end)) = last_TRF.Phase(end);


%==============================
function TRF = local_read(FN);
DD = fullfile(ZWOAEdatadir, 'calibration');
FN = fullfilename(FN, DD, '.probeloss');
if ~exist(FN,'file'),
    error(['File ''' FN ''' not found.']);
end
load(FN, '-mat'); % this loads varribale TRF






