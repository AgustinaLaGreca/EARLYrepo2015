function [DL, Dphi, Dt, binCompDt]=calibrate(EXP, Fsam, DAchan, Fstim, GivePhaseFac);
% calibrate - amplitude and phase transfer for acoustically flat transfer
%   [DL, Dphi, Dt, binCompDt]=calibrate(EXP, Fsam, chanStr, Fstim) computes 
%   the level term DL and phas e term Phi that compensate the acoustic
%   transfer at stimulus frequency Fstim determined during calibration, and
%   the group delay Dt of the DAC-to-acoustical path and the waveform delay
%   binCompDt that must be used to compensate interaural difference in Dt
%   values. The inputs are
%        EXP: experiment definition including the calibration data
%       Fsam: sample rate Hz
%    chanStr: DAchannel, L|R
%      Fstim: stimulus freq [Hz]. If none of the stim freqs are positive, 
%             and at least one stim freq is negative, the absolute value of
%             stimfreq is used, and anything outside the calibration range 
%             gets DL=0 dB and Dphi=0 cycle. The default behavior (i.e. 
%             when specifying all non-negative Fstim) is to throw an
%             error when any component of Fstim is outside the calibration
%             range of EXP.
%   The outputs are
%         DL: Level difference [dB] between numerical and acoustical
%             intensities. If W is the numerical sinusoid to be sent to the
%             DA converter, and the intensity of W is Lw dB_num,  then the
%             corresponding tone has a level of SPL = Lw-DL dB SPL. Here
%             dB_num is the power in dB of array W re RMS=1. Thus when
%             computing a tone at frequency Fstim, the amplitude of the
%             tone must be A = sqrt(2)*dB2A(X+DL) in order to realize a
%             sound pressure level of X dB SPL. The sqrt(2) factor comes
%             from the RMS of a unity-amplitude sinusoid, which is sqrt(1/2).
%       Dphi: Phase difference (cycle) between numerical and acoustical
%             waveforms. In order to get, at time t, a phase Phi_acoust of 
%             the acoustical waveform, the phase Phi of the numerical 
%             waveform W must equal Phi = Phi_acoust+Dphi. Importantly,
%             this phase difference does *not* include the overall group
%             delay Dt evaluated over the complete range calibration
%             frequencies. This is because one should not compensate a 
%             true time delay by phases of the components in cases in which
%             the envelope is constructed in the time domain (gating), as
%             this would create a discrepancy between the timing of the
%             waveform fine structure and that of the imposed gating.
%         Dt: global group delay [ms] of the calibration data. Estimate of 
%             the electrical+acoustical delay. If a wideband waveform is 
%             converted at time zero, the mean arrival time of its 
%             acoustical energy is at Dt. This group delay should not be
%             "compensated" when computing stimulus waveforms (see above 
%             comments on Dphi). The Dt parameter, however, is important 
%             at analysis time, when interpreting spike data and their 
%             latencies and group delays re the acoustical waveform.
%  binCompDt: zero if this DAchannel is lagging the other in terms of Dt;
%             minus the difference if this DAchannel is the leading one.
%             Delaying the waveform of this DAchannel by binDt will 
%             compensate the interaural asymmetry in latency caused by 
%             acoustical or electrical differences in the D/A-to-eardrum 
%             path. When only one channel is active (monaural calibration 
%             data), binCompDt is zero by convention.
%
%   [PhaseFac, dum, Dt, binCompDt]=calibrate(EXP, Fsam, chanStr, Fstim, 1)
%   returns the calibration as complex phase factors, i.e.,
%        PhaseFac = db2a(DL).*exp(2*pi*i*Dphi)
%
%   Note that both Dt and binDt are independent of Fstim, but may vary
%   with Fsam.
%          
%   When stimFreq is an array, both L and Phi are arrays of the same size.
%     
%   See also Experiment/edit, Experiment/interauralWBdelay. 

GivePhaseFac = arginDefaults('GivePhaseFac',0);
BeFreqTolerant = all(Fstim<=0) && any(Fstim<0); % see help text
Fstim = abs(Fstim);

% strip channel spec
DAchan= strrep(DAchan, '=Ipsi', ''); 
DAchan= strrep(DAchan, '=Contra', '');
[DAchan, Mess] = keywordMatch(DAchan,{'Left' 'Right'}, 'DAchan');
error(Mess);
ichan = strcmp(DAchan, {'Left' 'Right'});
ichanAvailable = strcmp(EXP.AudioChannelsUsed, {'Left' 'Right' 'Both'});
if isequal(3,ichanAvailable), ichanAvailable=[1 2]; end;
if ~ismember(ichan, ichanAvailable),
    error(['Audio channel ''' DAchan ''' not used in Experiment ''' name(EXP) '''.']);
end
Dt = nan;
switch EXP.CalibrationType,
    case 'Flat',
        SPL_1_Volt = str2num(strtok(EXP.Probes{ichan})); % dB SPL when RMS is 1 Volt
        DL = 0*Fstim - SPL_1_Volt;
        Dphi = 0*Fstim;
        [Dt, binCompDt] = deal([0 0],0);
    case 'Probe', % apply the cavity TRF fcn to the requested frequencies
        ProbeName = EXP.Probes{ichan};
        % get cached version, if it exists
        Cname = [mfilename '_probe_' ProbeName];
        CAV = MyFlag(Cname);
        if isempty(CAV), % no cached version yet. 
            CAV = cavityTrf(load(probetubecalib, EXP.Probes{ichan}));
            CAV = sparsify(CAV,0.001);
            MyFlag(Cname, CAV);
        end
        Ztr = eval(CAV, Fstim, BeFreqTolerant);
        DL = -A2dB(abs(Ztr));
        Dphi = -cangle(Ztr);
    case 'In situ',
        [idx, datenum]=calib_index(EXP);
        if isnan(datenum),
            error(['No ear calibration data found for experiment ' upper(name(EXP)) '.']);
        end
        Ename = [mfilename '_ear_' num2str(round(datenum*1e5))];
        EAR = MyFlag(Ename);
        if isempty(EAR), % no cached version yet. 
            EAR = transfer(load(earcalib, EXP, inf)); % most recent
            EAR = sparsify(EAR,0.001);
            MyFlag(Ename, EAR);
        end
        Ztr = eval(EAR(ichan), Fstim, BeFreqTolerant);
        DL = -A2dB(abs(Ztr));
        Dphi = -cangle(Ztr);
        Dt = getWBdelay(EAR);
    otherwise,
        error(['Unknown calibration Type ''' EXP.CalibrationType '''.']);
end
if isnan(Dt), % still to do
    Dt = calibWBdelay(EXP);
end
binCompDt = max(Dt)-Dt(ichan);
Dt = Dt(ichan);

if GivePhaseFac,
    DL = dB2A(DL).*exp(2*pi*i*Dphi);
    Dphi = [];
end

return








tau = [1.132 1.087]; % ms  [Left Right] 
ichan = strmatch(DAchan, {'Left', 'Right'});
Dt = tau(ichan);

% just cook something up. 1 dBnum ~ 1 Volt RMS ~ 90 dB
Fhi = 5e3*[1 pi^0.1]; 
DL = 12*exp(-((Fstim-1111)/820).^2)+17*exp(-((Fstim-3354)/2200).^2)+11*exp(-Fstim/Fhi(ichan));
DL = -(90+DL)+22;
if isempty(DL), Dphi=[];
else,
    Dphi = diff(db2a(DL)); Dphi = Dphi-mean(Dphi); Dphi=0.11*Dphi/max(abs(Dphi));
    if isempty(Dphi), Dphi=0; else, Dphi = Dphi([1 1:end]); end
end

if isequal('Both', EXP.AudioChannelsUsed),
    % binCompDt is only nonzero if ichan is the faster one
    binCompDt = tau(3-ichan)-tau(ichan); % delay of the opposite chan re delay of this one this one
    binCompDt = max(0,binCompDt); % if binCompDt<0, the other chan is faster -> set binCompDt to zero.
else, binCompDt = 0; % see help
end

%================
Dt = 0*Dt;
binCompDt = 0;
Dphi = 0*Dphi;




