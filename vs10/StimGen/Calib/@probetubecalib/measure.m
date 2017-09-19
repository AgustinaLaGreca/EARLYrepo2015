function [PTC, Mess, ErrField] = measure(PTC, P, figh);
% Probetubecalib/measure - measure transfer function of probe tube
%    S = ProbeCalib(Probetubecalib(), Param, H) performs a system calibration
%    for measuring the acoustic transfer function of a probe tube. A cavity
%    is used with a microphone in place of the eardrum. The outputs of the microphone
%    amplifier are connected to the analog inputs of the RX6 as follows
%      Ch 1: probe microphone (typically a 1/2" B&K)
%      Ch 2: cavity mic (the one mimicking the eardum, say a 1/4" B&K)
%    Param is the struct obtained from calling GUIval on the GUI produced
%    by ProbeTubeCalibration. The optional input argument H is a figure
%    handle of a GUI. If it is defined, GUImessage is used to report the
%    progress of the recording. Example:
%              PeakAmp: 1
%          PeakAmpUnit: 'V'
%                Atten: 37.8000
%            AttenUnit: 'dB'
%               DAchan: '1'
%           DAchanUnit: ''
%                 Fmin: 50
%             FminUnit: 'Hz'
%                 Fmax: 45000
%             FmaxUnit: 'Hz'
%                Speed: 1
%            SpeedUnit: 'octave/s'
%               handle: [1x1 struct]
%              GUIgrab: [1x1 struct]
%   If something is wrong with any of these parameters, an error is thrown.
%   T is a struct with fields 
%         Probe: transfer of D/A to probe mic
%        Cavity: transfer of D/A to cavity mic
%          Tube: transfer from probe to cavity (that is inverse(Probe)*Cavity)
%  
%   [T, Mess, ErrField]=ProbeCalib(Param) suppresses any errors with the
%   Parameters and returns a void T and message Mess instead. ErrField is 
%   the name of the offending field in Param. This can be used to highlight 
%   the appropriate field in the GUI. 
%
%   NOTE 
%   The first input arg is a dummy for Methodizing this function.
%
%   See also ProbeTubeCalibration, Methodizing.

MaxAbsDA = 9.9; % Volt maxabs D/A
if nargin<3, figh=[]; end % no reporting to a GUI
MaxNsamSweep = 2e5; % indiv sweeps should not last much more than MaxNsamSweep samples

doReport = isSingleHandle(figh);
if isequal('stop',P),
    local_measure('stop');
    return;
else, % reset stop flag in local_measure
    local_measure('nostop');
end

SupressError = nargout>1;
Fsam = 1e3*RX6sampleRate(130); % Hz

% Check the parameter values
Noct = log2(P.Fmax/P.Fmin);
Dur = abs(Noct/P.Speed); % total net stim duration in s
if P.PeakAmp>10,
    Mess = 'PeakAmp may not exceed 10 V for RX6 D/A converters.';
    ErrField = 'PeakAmp';
elseif P.Atten>120,
    Mess = 'Analog attenuation may not exceed 120 dB for PA5s.';
    ErrField = 'Atten';
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

if doReport, GUImessage(figh, 'Preparing calib stimulus...'); end
sys3PA5( P.Atten, nan); %changed by EVE 13/08/2015 old version = sys3PA5(P.Atten);

T0 = transfer('voltage', 'pressure', '0 volt', '20 uP', 'dBV', 'dB SPL');

Nsam = Dur*Fsam;
Nband = ceil(Nsam/MaxNsamSweep); % indiv sweeps should not last much more than MaxNsamSweep samples
Flim = logispace(P.Fmin, P.Fmax, Nband+1);
Overlap = min(sqrt(2), sqrt(Flim(2)/Flim(1))); % overlap of adjacent bands in octaves
Fstart = round([P.Fmin, Flim(2:end-1)/Overlap]); % anticipate 1/2-octave overlap of consecutive sweeps
Fend = round([Flim(2:end-1)*Overlap P.Fmax]); % 1/2-octave overlap of consecutive sweeps

%----Probe  (chan 1): 1/2" B&K---------
HW1 = struct('DAdevice', 'RX6', 'ADdevice', 'RX6', 'DAchan', str2num(P.DAchan), 'ADchan', 1);
Sens1 = struct('Stim', dB2A(P.Atten), 'Resp', P.ProbeSens/dB2A(94)); 
Amp1 = struct('Freq', [1 1e6], 'Amp_stim', P.PeakAmp, 'SNR', inf, 'MaxAbsDA', MaxAbsDA);
Descr1 = ['Calib [Volt DAC -> dB SPL] via chan 1 of B&K mic amplifier (' num2str(P.ProbeSens) ' V/Pa)'];
for iband=1:Nband,
    Freq1(iband) = struct('Fsam', Fsam, 'Fmin', Fstart(iband), 'Fmax', Fend(iband), 'SmoothBand', 0.005);
end
Tim1 = struct('Speed', P.Speed, 'Pdur', 200);

%----chan 2: 1/4" B&K---------
HW2 = struct('DAdevice', 'RX6', 'ADdevice', 'RX6', 'DAchan', str2num(P.DAchan), 'ADchan', 2);
Sens2 = struct('Stim', dB2A(P.Atten), 'Resp', P.CavitySens/dB2A(94)); 
Amp2 = Amp1;
Descr2 = ['Calib [Volt DAC -> dB SPL] via chan 2 of B&K mic amplifier (' num2str(P.CavitySens) ' V/Pa)'];
Freq2 = Freq1;
Tim2 = Tim1;


%  ========== the measurements ==============
Probe = local_measure(T0, HW1, Freq1, Tim1, Sens1, Amp1, {Descr1}, figh, 'Probe: recording');
if isvoid(Probe), return; end % interrupted
Cavity = local_measure(T0, HW2, Freq2, Tim2, Sens2, Amp2, {Descr2}, figh, 'Cavity: recording');
if isvoid(Cavity), return; end % interrupted
Tube = inverse(Probe)*Cavity; % note the order, yielding acoustic -> acoustic transfer
Tube = setWBdelay(Tube,'wflat');
Tube = description(Tube, 'Inverse probe tube TRF from 1/2" mike to eardrum');
PTC = probetubecalib(Probe, Cavity, Tube); % collect T into Probetubecalib object


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
        sys3PA5(120);
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
    T(1) = patch(T(1),T(iband)); % 1+2, (1+2)+3 ...
end
T = setWBdelay(T(1), 'wflat');



function local_hold;
hh = findobj('type','axes'); 
set(gcf,'NextPlot', 'add'); 
set(hh,'NextPlot', 'add');







