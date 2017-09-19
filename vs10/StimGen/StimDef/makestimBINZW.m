function P2=makestimBINZW(P);
% MakestimRF - stimulus generator for BINZW stimGUI (bin tone complex)
%    P=MakestimBINZW(P), where P is returned by GUIval, generates the stimulus
%    specified in P. MakestimZW is typically called by StimGuiAction when
%    the user pushes the Check, Play or PlayRec button.
%    MakestimRF does the following:
%        * Complete check of the stimulus parameters and their mutual
%          consistency, while reporting any errors
%        * Compute the stimulus waveforms
%        * Computation and broadcasting info about # conditions, total
%          stimulus duration, Max SPL, etc.
%
%    MakestimBINZW renders P ready for D/A conversion by adding the following 
%    fields to P
%           Fsam: sample rate [Hz] of all waveforms. This value is
%                 determined by carrier & modulation freqs, but also by
%                 the Experiment definition P.Experiment, which may 
%                 prescribe a minimum sample rate needed for ADC.
%           Fcmp: frequencies [Hz] of all the components in a column array
%            ITD: ITDs [ms] of all the conditions
%       Waveform: Waveform object array containing the samples in SeqPlay
%                 format.
%    Attenuation: scaling factors and analog attuater settings for D/A
%   Presentation: struct containing detailed info on stimulus order,
%                 broadcasting of D/A progress, etc.
% 
%   See also Waveform/maxSPL, Waveform/play, sortConditions, 

% Note: fixed ITD  sign bug, July 31, 2012. Experiments RG12448 and before
% have inverted ITDs!

P2 = []; % a premature return will result in []
if isempty(P), return; end
figh = P.handle.GUIfig;
EXP = P.Experiment;

% check & convert params. Note that helpers like maxSPL
% report any problems to the GUI and return [] or false in case of problems.

%     i_________ToneComplex_______: '_____ToneComplex_____'
%                          LowFreq: 100
%                      LowFreqUnit: 'Hz'
%                         HighFreq: 2000
%                     HighFreqUnit: 'Hz'
%                            Ncomp: 30
%                        NcompUnit: ''
%                          FreqTol: 10
%                      FreqTolUnit: '%'
%                           ZWseed: 697192199
%                       ZWseedUnit: ''
somethingwrong=1;
if P.LowFreq<EXP.minStimFreq,
    GUImessage(figh, {'Low frequency violates min stim frequency'...
        ['of ' num2str(EXP.minStimFreq) ' Hz']},'error', 'LowFreq');
elseif P.HighFreq>EXP.maxStimFreq,
    GUImessage(figh, {'High frequency exceeds max stim frequency'...
        ['of ' num2str(EXP.minStimFreq) ' Hz']},'error', 'HighFreq');
elseif P.LowFreq>P.HighFreq,
    GUImessage(figh, {'Low frequency exceeds High frequency'}, 'error', {'LowFreq' 'HighFreq'});
elseif P.Ncomp<2,
    GUImessage(figh, {'Too few (<2) freq components'}, 'error', {'Ncomp'});
elseif P.Ncomp>100,
    GUImessage(figh, {'Too many (>100) freq components'}, 'error', {'Ncomp'});
elseif P.FreqTol>50,
    GUImessage(figh, {'Frequency tolerance exceeds max of 50%.'}, 'error', {'FreqTol'});
else, % passed the tests..
    somethingwrong=0;
end
if somethingwrong, return; end

% freq range & sample rate
F_reg = linspace(P.LowFreq, P.HighFreq, P.Ncomp).'; % regular spacing
AvgDF = diff(F_reg(1:2)); % average freq spacing of components
Freq_AbsTol = 0.01*P.FreqTol*AvgDF; % absolute value of tolerance (one-sided)
P.Fsam = sampleRate(P.HighFreq+Freq_AbsTol, EXP); % sample rate in Hz
ZwM = find(messenger(), figh, 'ZwuisMess');
report(ZwM, ['Average spacing: ' num2str(round(AvgDF)) ' Hz']);
% reportFsam(figh, P.Fsam/1e3); % in kHz
%             i_________ITD_______: '_____ITD_____'
%                         startITD: -1
%                          stepITD: 0.1000
%                           endITD: 1
%                        AdjustITD: 'step'
[P.ITD, Mess] = EvalStepper(P.startITD, P.stepITD, P.endITD, 'linear');
if isequal('nofit', Mess),
    GUImessage(figh, {'No integer # ITD steps fit in bounds. Change tolerance or ITD step specs.'}, 'error', ...
        {'ITDstart' 'ITDstep' 'ITDend' });
elseif isequal('largstep', Mess),
    GUImessage(figh, {'Stepsize too big to fit in bounds.'}, 'error', ...
        {'ITDstart' 'ITDstep' 'ITDend' });
end
if ~isempty(Mess), return; end
Ncond = numel(P.ITD);
if Ncond>100,
    GUImessage(figh, {'Too many (>100) ITD values. Increase step or decrease range.'}, 'error', ...
        {'ITDstart' 'ITDstep' 'ITDend' });
    return;
end
%             i_________Dur_______: '_____Dur_____'
%                         BurstDur: 500
%                          RiseDur: 10
%                          FallDur: 10
P.ITDtype = 'ongoing'; % 
if ~EvalDurPanel(figh, P, Ncond), 
    return;
end
%            i_________Pres_______: '_____Pres_____'
%                              ISI: 100
%                             Nrep: 20
%                         Baseline: 1000
%                         Grouping: 'by condition'
%                            Order: 'Forward'
%                            RSeed: 402202720

% Compute the waveforms (P.Waveform) 
P = local_waveform(P, F_reg, Freq_AbsTol);

% report duration
totBaseline = sum(SameSize(P.Baseline,[1 1])); % sum of pre- & post-stim baselines
Ttotal=ReportPlayTime(figh, Ncond, P.Nrep, P.ISI, totBaseline);


%       i_________Levels_______: '_____Levels_____'
%                           SPL: 60
%                           DAC: 'Left=Ipsi'
P = sortConditions(P, 'ITD', 'ITD' , 'ms', 'Lin');

% Levels and active channels (must be called *after* adding the baseline waveforms)
[mxSPL P.Attenuation] = maxSPL(P.Waveform, P.Experiment);
okay = CheckSPL(figh, P.SPL, mxSPL, P.ITD, '', {'SPL'});
if ~okay, return; end

P2=P;

%===================================================
function P = local_waveform(P, F_reg, Freq_AbsTol);
% F_reg is regular frequency array. Freq_Abstol is max deviation from F_reg
SetRandState(P.ZWseed); 
Fbound = [F_reg-Freq_AbsTol, F_reg+Freq_AbsTol];
w = rand(P.Ncomp,1); 
P.Fcmp = sum(Fbound.*[w 1-w],2);
dt = 1e3/P.Fsam;
Ncmp = numel(P.Fcmp);
Ncond = numel(P.ITD);
P.StartPhase = rand(1,Ncmp); % starting phases of zwuis components
AllChan = channelSelect(P.DAC, 'LR');
chanDelay = ITD2delay(P.ITD, P.Experiment); % ITD-realizing delays for indiv D/A channels
for DAchan=AllChan,
    ichan = 1+double(P.DAC(1)=='B' && DAchan=='R');
    for icond = 1:Ncond,
        [chdel, SPL, risedur, falldur] = channelSelect(DAchan, chanDelay(icond,:), P.SPL, P.RiseDur, P.FallDur);
        ITDphaseShift = -1e-3*chdel*P.Fcmp; % delay -> phase in cycles 
        [DL, Dphi] = calibrate(P.Experiment, P.Fsam, DAchan, P.Fcmp);
        LinAmp = dB2A(SPL+DL)*sqrt(2); % calibrated linear amplitude
        ph0 = P.StartPhase(:)+Dphi(:)+ITDphaseShift; % calibrated starting phase, including ITD
        % function T = tonecomplex(Amp, Freq, Ph0, Fsam, Dur);
        w = tonecomplex(LinAmp, P.Fcmp, ph0, P.Fsam, P.BurstDur+dt); % one extra sample to avoid compulsive complaints by exactGate
        %function X = ExactGate(X, Fsam, BurstDur, Tstart, RiseDur, FallDur);
        w = ExactGate(w, P.Fsam, P.BurstDur, 0, P.RiseDur, P.FallDur);
        Param = CollectInStruct(LinAmp, ph0);
        P.Waveform(icond, ichan) = Waveform(P.Fsam, DAchan, max(abs(w)), P.SPL, Param, {w}, 1);
        % trailing zeros
        P.Waveform(icond, ichan) = AppendSilence(P.Waveform(icond, ichan), P.ISI);
    end
end

function  [Rise, Fall] = local_ramps(w, Nsam);
Wrise = sin(linspace(0,pi/2,Nsam).').^2;
Wfall = flipud(Wrise);
Rise = Wrise.*w(end+1+(-Nsam:-1));
Fall = Wfall.*w(end+1+(-Nsam:-1));
 





