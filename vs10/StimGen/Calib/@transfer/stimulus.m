function [Zstim, InstFreq, StimInfo] = stimulus(T, flag);
% transfer/stimulus - complex stimulus waveform for calibration measurement
%    [Zstim, InstFreq, StimInfo] = stimulus(T) returns the complex waveform 
%    Zstim (unscaled for D/A sensitivity), instantaneous frequency InstFreq
%    in Hz, and various stimulus parameters in struct StimInfo.
%
%    stimulus(T, 'trim') restricts the stimulus to the sweeping part.
%
%    See Transfer, Transfer/measure.

if nargin<2, flag=''; end

CP = T.CalibParam;
if isempty(CP),
    error('Calibration params of Transfer object T not specified.');
end

dt = 1e3/CP.Freq.Fsam; % ms sample period
Noct = log2(CP.Freq.Fmax/CP.Freq.Fmin);
SweepDur = 1e3*Noct/abs(CP.Tim.Speed); % sweep duration in ms
NsamSweep = 10*round(SweepDur/dt/10); % # samples of sweep
NsamPlateau = 10*round(CP.Tim.Pdur/dt/10); % # samples of plateaus
NsamRamp = round(NsamPlateau/2); % # samples of cos^2 gating
InstFreq = logispace(CP.Freq.Fmin, CP.Freq.Fmax, NsamSweep).'; % col vector holding instantaneous freq of sweep in Hz
InstFreq = [ones(NsamPlateau,1)*CP.Freq.Fmin; InstFreq; ones(NsamPlateau,1)*CP.Freq.Fmax]; % provide plateaus
if numel(CP.Amp.Amp_stim)==1, % constant output
    InstAmp = CP.Amp.Amp_stim;
else, % stimulus amplitude varies with frequency; use interpolation to get instantaneous values
    InstAmp = interp1(CP.Amp.Freq, CP.Amp.Amp_stim, InstFreq);
end
if any(isnan(InstAmp)),
    error('Error in specification of frequency range/amplitude of TRF stimulus.');
end
InstPhase = cumsum((2*pi*1e-3*dt)*InstFreq); % instantaneous phase in radians
Zstim = InstAmp.*exp(i*InstPhase); % complex stimulus
Zstim = simplegate(Zstim, NsamRamp); % provide cos^2 ramps
if CP.Tim.Speed<0, % downward sweep
    Zstim = conj(flipud(Zstim));
    InstFreq = flipud(InstFreq);
end
StimInfo = CollectInStruct(NsamSweep, NsamPlateau, NsamRamp);
if isequal('trim', flag),
    isweep = StimInfo.NsamPlateau+(1:StimInfo.NsamSweep); % indices of sweeping part of stim
    [Zstim, InstFreq] = deal(Zstim(isweep), InstFreq(isweep));
end





