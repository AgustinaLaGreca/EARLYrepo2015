function [Trf, Mess, ErrField] = measure(EC, ichanDA, ichanAD, P, figh);
% Earcalib/measure - measure in situ transfer function of ear canal
%    Trf = measure(EC, ichanDA, ichanAD, Param, H) performs an ear calibration
%    for measuring the acoustic transfer function Trf from the specified
%    channel (# ichanDA) of the D/A device to the ear canal. The microphone
%    signal is connect to channel # ichanAD. 
%    These data are typically used for calibration during an experiment.
%    The output of the probe microphone amplifier must be connected to 
%    the analog input #1 of the RX6. At the moment of the call, EC must
%    contain at least the Experiment and activeDAchan fields.
%
%    Param is the struct obtained from calling GUIval on the GUI produced
%    by earcalib/GUI. The optional input argument H is a figure
%    handle of that GUI. If it is defined, GUImessage is used to report the
%    progress of the recording. Example of Param:
%         SpeakerName_1: 'shure1'
%     SpeakerName_1Unit: ''
%         SpeakerName_2: 'shure2'
%     SpeakerName_2Unit: ''
%           ProbeName_1: 'needle_1'
%           ProbeSens_1: 1
%       ProbeSens_1Unit: 'V/Pa'
%           ProbeName_2: 'needle_2'
%           ProbeSens_2: 1
%       ProbeSens_2Unit: 'V/Pa'
%                  Fmin: 50
%              FminUnit: 'Hz'
%                  Fmax: 45000
%              FmaxUnit: 'Hz'
%                 Speed: 2
%             SpeedUnit: 'octave/s'
%                MaxSPL: 85
%            MaxSPLUnit: 'dB SPL'
%   If something is wrong with any of these parameters, an error is thrown.
%  
%   [T, Mess, ErrField]=EarCalib(...) suppresses any errors with the
%   Parameters and returns a void T and message Mess instead. ErrField is 
%   the name of the offending field in Param. This can be used to highlight 
%   the appropriate field in the GUI. 
%
%   See also earcalibration, Methodizing.
MaxAbsDA = 9.9; % Volt maxabs D/A
MaxNsamSweep = 2e5; % indiv sweeps should not last much more than MaxNsamSweep samples

figh = arginDefaults('figh', []); % no reporting to a GUI

doReport = isSingleHandle(figh);
if isequal('stop',ichanDA),
    local_measure('stop');
    return;
else, % reset stop flag in local_measure
    local_measure('nostop');
end

Trf = transfer(); % void default return arg
SupressError = nargout>1;
Fsam = 1e3*RX6sampleRate(130); % Hz
chanDAstr = num2str(ichanDA);

% get probetube Trf
PTC_FN = P.(['ProbeName_' chanDAstr]); % filename of probe tube calib
% Check the parameter values
Noct = log2(P.Fmax/P.Fmin);
Dur = abs(Noct/P.Speed); % total net stim duration in s
if ~exist(fullfilename(probetubecalib(), PTC_FN), 'file'),
    Mess = ['Cannot file probe tube file ''' PTC_FN  '''.'];
    ErrField = ['ProbeName_' chanDAstr];
elseif P.Fmin<25,
    Mess = 'Fmin must be >= 25 Hz.';
    ErrField = 'Fmin';
elseif P.Fmax<25,
    Mess = 'Fmax must be >= 25 Hz.';
    ErrField = 'Fmax';
elseif P.Fmin>50e3,
    Mess = 'Fmin must be <= 50 kHz.';
    ErrField = 'Fmin';
elseif P.Fmax>50e3,
    Mess = 'Fmax must be <= 50 kHz.';
    ErrField = 'Fmax';
elseif P.Fmin>=P.Fmax,
    Mess = 'Fmin may not be larger than Fmax.';
    ErrField = {'Fmin' 'Fmax'};
elseif abs(P.Speed)>5,
    Mess = 'Absolute value of Speed may not exceed 5 octave/s.';
    ErrField = 'Speed';
elseif Dur>30,
    Mess = 'Calibration takes too long (>30 s). Increase Speed';
    ErrField = 'Speed';
else, % okay
    Mess = '';
    ErrField = '';
end
% Error handling depends on # out args. See help text.
if ~isempty(Mess),
    if SupressError, return;
    else, error(Mess);
    end
end

% convert max SPL into max DA ampl using the probe
PTC = load(probetubecalib, PTC_FN); 
CavityTransfer = cavityTrf(PTC); % this is TRF (DA dBV -> cavity dBSPL)
TubeTransfer = tubeTrf(PTC); % this is acoust TRF (probe -> ear canal)
MaxSPL_at_Probe = P.MaxSPL - magnitude(TubeTransfer); % max SPL @ probe
% estimated DA peak voltage to reach max dB SPL in ear canal (~cavity)
guessed_MaxDA_amp = dB2A(P.MaxSPL-magnitude(CavityTransfer))*sqrt(2);
DA_amp = min(guessed_MaxDA_amp, 0.9*MaxAbsDA); % safety margin of 0.9 to stay away from clipping the D/A
freq = frequency(TubeTransfer);
% subplot(2,1,1); plot(freq, MaxSPL_at_Probe);
% subplot(2,1,2); plot(freq, guessed_MaxDA_amp); ylog125
% xplot(freq, DA_amp, 'r');
if ~betwixt(1.000000001*P.Fmin, min(freq), max(freq)),
    Mess = 'Min Frequency lies outside the range of the probe tube calibration data.';
    ErrField = 'Fmin';
elseif ~betwixt(P.Fmax/1.000000001, min(freq), max(freq)),
    Mess = 'Max Frequency lies outside the range of the probe tube calibration data.';
    ErrField = 'Fmax';
end
if ~isempty(Mess),
    if SupressError, return;
    else, error(Mess);
    end
end

% prepare recording
if doReport, GUImessage(figh, 'Preparing calib stimulus...'); end
T0 = transfer('voltage', 'pressure', '0 volt', '20 uP', 'dBV', 'dB SPL');
% divide recordings over frequency bands
Nsam = Dur*Fsam;
Nband = ceil(Nsam/MaxNsamSweep); % indiv sweeps should not last much more than MaxNsamSweep samples
Flim = logispace(P.Fmin, P.Fmax, Nband+1);
Overlap = min(sqrt(2), sqrt(Flim(2)/Flim(1))); % overlap of adjacent bands in octaves
Fstart = round([P.Fmin, Flim(2:end-1)/Overlap]); % anticipate 1/2-octave overlap of consecutive sweeps
Fend = round([Flim(2:end-1)*Overlap P.Fmax]); % 1/2-octave overlap of consecutive sweeps

probeSens = P.(['ProbeSens_' chanDAstr]);
HW = struct('DAdevice', 'RX6', 'ADdevice', 'RX6', 'DAchan', ichanDA, 'ADchan', ichanAD);
Sens = struct('Stim', 1, 'Resp', probeSens/dB2A(94)); 
Amp = struct('Freq', freq, 'Amp_stim', DA_amp, 'SNR', inf, 'MaxAbsDA', MaxAbsDA);
Descr = ['Calib [Volt DAC ' chanDAstr  ' -> dB SPL] via chan ' num2str(ichanAD) ' of B&K mic amplifier (' num2str(probeSens) ' V/Pa)'];
for iband=1:Nband,
    Freq(iband) = struct('Fsam', Fsam, 'Fmin', Fstart(iband), 'Fmax', Fend(iband), 'SmoothBand', 0.005);
end
Tim = struct('Speed', P.Speed, 'Pdur', 200);



%  ========== the measurements ==============
if AttenPresent, sys3PA5(0); end
Trf = local_measure(T0, HW, Freq, Tim, Sens, Amp, {Descr}, figh, 'recording');
if isvoid(Trf), return; end % interrupted


%============================================
function T=local_measure(T0, HW, Freq, Tim, Sens, Amp, Descr, figh, ADname);
%dsize(HW, Freq, Tim, Sens, Amp, Descr)
persistent doStop
if isempty(doStop), doStop=0; end; 
if isequal('stop',T0),
    doStop = 1;
    return;
elseif isequal('nostop',T0),
    doStop = 0;
    return;
end
doReport = isSingleHandle(figh);
[HW, Freq, Tim, Sens, Amp, Descr] = SameSize(HW, Freq, Tim, Sens, Amp, Descr);
Nband = numel(Freq);
for iband = 1:Nband,
    if doReport,
        GUImessage(figh,[ADname ' band ' num2str(iband) '/' num2str(Nband) '.'])
    end
    pause(0.05); % allow interrupts to interfere
    if doStop,
        if AttenPresent, sys3PA5(120); end
        seqplayhalt;
        T = transfer; % void default return value to signal interrupt
        return;
    end
    T(iband) = measure(T0, HW(iband), Freq(iband), ...
        Tim(iband), Sens(iband), Amp(iband), Descr(iband));
end
T = setWBdelay(T, max(getWBdelay(T)));
% patch all bands into T(1)
for iband = 2:Nband,
    T(1) = patch(T(1),T(iband)); % 1+2 -> (1+2)+3 ...
end
T = setWBdelay(T(1), 'wflat');



function local_hold;
hh = findobj('type','axes'); 
set(gcf,'NextPlot', 'add'); 
set(hh,'NextPlot', 'add');







