function [T, Rec] = measure(T, Dev, Freq, Tim, Sens, Amp, Descr);
% transfer/measure - measure transfer function
%    T = measure(T, H, F, T, S, A, Descr) measures a transfer function 
%    using a log sweep stimulus and stores it in a previously defined Transfer 
%    object T.
%    Inputs:
%      Dev: hardware info. Struct with fields
%           DAdevice:  device name for D/A. Currently only RX6.
%           DAdevice:  device name for A/D. Currently only RX6.
%             DAchan:  D/A channel. Default is channel 1 (Analog-OUT-1).
%             ADchan:  A/D channel. Default is channel 1 (Analog-IN-1).
%      F: frequency info. Struct with fields
%               Fsam:  sample rate in Hz
%               Fmin:  lowest freqency in Hz
%               Fmax:  highest freqency in Hz
%      T: timing info. Struct with fields
%             Speed: FM speed in Octaves/s. Negative sweeps result in
%                    downward sweeps.
%              Pdur: duration (ms) of the plateaus before & after the sweep.
%                    cos^2 ramps are half this duration.
%      S: sensitivity info. Struct with fields
%              Stim: numerical conversion factor relating the RMS of the
%                    waveform sent to the D/A device to the reference in
%                    T.Ref_stim. When the RX6 is used for D/A without an 
%                    attenuator, this factor equals 1 V/V. If an analog 
%                    attenuation of, say, 40 dB is used, the conversion
%                    factor is 100 V/V.
%                    As a mnemonic, think of the conversion factor as 
%                      "Sens.Stim = DA_amplitude / stimulus amplitude"
%              Resp: numerical conversion factor relating the RMS of the
%                    response (waveform obtained from the A/D) to the resp 
%                    reference Ref_resp of T. For instance, for a microphone
%                    signal, when the sensitivity is 200 mV/Pa, and the
%                    recorded waveform is obtained from the A/D of the RX6,
%                    then the conversion factor equals 0.200/a2db(94) V/Pa, 
%                    because 1 Pa corresponds to 94 dB SPL and the signal 
%                    obtained from the RX6 is given in volts. 
%                    As a mnemonic, think of the conversion factor as 
%                      "Sens.Resp = AD_amplitude / response amplitude"
%       A: amplitude info. Struct with fields
%              Freq: frequencies [Hz] to which the amplitudes (see next field)
%                    correspond. [] for constant amplitude.
%          Amp_stim: amplitude [dB] of the stimulus expressed re Ref_stim. 
%                    Either a constant (flat spectrum) or an array 
%                    corresponding to the frequencies in A.Freq, from which
%                    the amplitude profile of the calibration stimulus will 
%                    be obtained by linear interpolation.
%          MaxAbsDA: maximum absolute value of D/A output, e.g. 10 for RX6.
%       D: Description of the calibration in char string. 

% test code @bottom

persistent lastTC

if isvoid(T),
    error('Cannot use a void Transfer object to store a transfer function.');
end

if length(Dev.DAdevice)<3 || ~isequal('RX6', upper(Dev.DAdevice(1:3))),
    error('Transfer/measure is currently implemented only for the RX6 as D/A device.');
end
if length(Dev.ADdevice)<3 || ~isequal('RX6', upper(Dev.ADdevice(1:3))),
    error('Transfer/measure is currently implemented only for the RX6 as A/D device.');
end
if ~isequal(Dev.DAdevice, Dev.ADdevice),
    error('D/A and A/D device must be the same RX6.');
end
if ~ismember(Dev.ADchan,[1 2]),
    error('Analog input channel # must be 1 or 2; forget TDT bullshit about channels 128 & 129.');
end

TmaxLag = 50; % ms max shift for overall delay between DA and AD

% get true sample rate by initializing the Seqplay circuit
[SP, RCLD] = seqplayinit(Freq.Fsam/1e3, Dev.DAdevice, '#Default', 1); % default seqplay circuit; 1=recycle circuit if possible
Fsam = 1e3*SP.Fsam; % kHz-> Hz
dt = 1e3/Fsam; % ms sample period

% perform temporal calibration DA->DA, i.e., determine delay between DA and
% DA when directly connected.
AnaAtten = sys3PA5get; % store attenuator settings. Note that sys3PA5 returns [] if no attenuators are connected.
if RCLD && ~isempty(lastTC),
    TC = lastTC;
else,
    TC = DA_AD_TimingCalibrateRX6(Dev.DAdevice, rechardware(experiment)); % this measures the DSP delay DAC->ADC (which is subject to time jitter)
    lastTC = TC;
end
% sys3PA5(AnaAtten);
sys3PA5(AnaAtten, nan); %changed by EVE 13/08/2015
drawnow;

% store stimulus params in T, so we can use Transfer/Stimulus method
Date = datestr(now);
Location = [Where '/' CompuName];
createdBy = [class(T) '/' mfilename];
MaxC = nan; % temp value (see below)
Freq.Fsam = Fsam;
T.CalibParam = CollectInStruct(Dev, Freq, Tim, Sens, Amp, AnaAtten, ...
    '-', TmaxLag, '-', TC, MaxC, '-', Date, Location, createdBy);

% prepare calibration stimulus
[Zstim, InstFreq, StimInfo] = stimulus(T);
Zstim = Zstim*Sens.Stim; % apply conversion factor stimulus -> D/A

%dplot(dt, real(Zstim))
%dplot(dt, InstFreq)

maxAbsSam = max(abs(Zstim(:)));
if maxAbsSam>Amp.MaxAbsDA,
    maxAbsSam
    error('Requested amplitude of calib stimulus too high to play without clipping.');
end

% ===========prepare D/A and A/D and give it a go
if isequal(1, Dev.DAchan),
    seqplayupload({real(Zstim)}, {});
    seqplaylist(1,1);
elseif isequal(2, Dev.DAchan),
    seqplayupload({}, {real(Zstim)});
    seqplaylist([],[], 1,1);
else,
    error('DA channel must be either 1 or 2.');
end
sys3setpar(1, 'DoGrab_ADC', Dev.ADdevice); % enable ADC
sys3trig(5,Dev.ADdevice);  % reset AD buffer counts
SPinfo = seqplayinfo;
NsamTot = SPinfo.NsamTotPlay;
seqplaygo;
% Mess = seqplayDAwait(1.5*NsamTot*dt); % timeout should be sufficient to play the whole stim
% error(Mess);

% ===========grab samples. We need to use double buffering, as the stimulus 
%             length easily exceeds the AD buffer size
ch_str = num2str(Dev.ADchan); 
DataTag = ['ADC_' ch_str];
IndexTag = ['Nsam_ADC_' ch_str];
Rec = sys3doubleBufRead(0, DataTag, IndexTag, numel(Zstim), Dev.ADdevice); % 0: no triggering by sys3doubleBufRead
Rec = (Rec - mean(Rec))/Sens.Resp; % reject DC; account for sensititivity factor
% filter out the sweep band; in passing make Rec a complex analytic signal
Rec = local_filt(Rec, Freq.Fsam, Freq.Fmin*0.9, Freq.Fmax*1.1);
% ============determine wideband delay and use it to time-align Stim & Rec
NsamMaxLag = round(TmaxLag/dt);
%[T.CalibParam.MaxC, iLag]=maxcorr(Rec, Zstim, NsamMaxLag, 'abs');
[T.CalibParam.MaxC, iLag]=maxcorr(real(Rec), real(Zstim), NsamMaxLag, 'abs');
T.CalibParam.StimInfo = StimInfo;
iLag = iLag-1; % make sure that residual delay is non-negative
Rec = circshift(Rec, -iLag);
%dplot(dt, real(Zstim)); xdplot(dt, real(Rec), 'r');
WBdelay = iLag*dt-TC.Lag_ms;  % correct for DA->AD delay of hardware

% ======get trf fun from stim & rec; smoothen a little
T.Ztrf = Sens.Stim*Rec./Zstim;
T.Ztrf = simplegate(T.Ztrf, StimInfo.NsamRamp);
T.Ztrf = local_LP(T.Ztrf, Amp.SNR);

doct = 1e-3*dt*abs(Tim.Speed); % spacing in octaves between neighboring points of TRF
doDec = (Freq.SmoothBand/doct)>100;
if doDec, 
    T.Ztrf = smoothen(T.Ztrf,20);
    T.Ztrf = 0.5*(T.Ztrf(5:10:end)+T.Ztrf(6:10:end));
    InstFreq = 0.5*(InstFreq(5:10:end)+InstFreq(6:10:end));
    isw0 = find(diff(InstFreq)>0,1,'first')-1;
    isw1 = find(diff(InstFreq)>0,1,'last');
    StimInfo.NsamSweep = isw1-isw0+1;
    StimInfo.NsamPlateau = isw0;
end
T.Ztrf = smoothen(T.Ztrf, Freq.SmoothBand, doct);

isweep = StimInfo.NsamPlateau+(1:StimInfo.NsamSweep); % indices of sweeping part of stim & rec
[T.Ztrf, T.Freq] = deal(T.Ztrf(isweep), InstFreq(isweep));
T.Description = Descr;
T.WB_delay_ms = WBdelay;



%====================================================
function X = local_filt(X, Fsam, Flo, Fhi);
% bandpass filtering; only retain pos freqs
%dsize(X)
X = fft(X);
Nsam = numel(X);
df = Fsam/Nsam;
%Flo, Fhi
ilo = 1+floor(Flo/df);
ihi = ceil(Fhi/df);
X = [zeros(ilo-1,1); X(ilo:ihi); zeros(Nsam-ihi,1)] ;
X = 2*ifft(X); % factor two compensates trashing neg freqs; now real(X) is equal to the original, real, X in the passband

function Trf = local_LP(Trf, SNR);
Trf = fft(Trf);
SP = A2dB(abs(Trf));
irej = find(SP<max(SP)-SNR, 1, 'first');
if ~isempty(irej),
    Trf(irej+1:end+1-irej) = 0;
end
Trf = ifft(Trf);



