function P2=makestimZWSPL(P);
% MakestimZWSPL - stimulus generator for ZWSPL stimGUI (bin tone complex)
%    P=MakestimZWSPL(P), where P is returned by GUIval, generates the stimulus
%    specified in P. MakestimZWSPL is typically called by StimCheck when
%    the user pushes the Check, Play or PlayRec button.
%    MakestimZWSPL does the following:
%        * Complete check of the stimulus parameters and their mutual
%          consistency, while reporting any errors
%        * Compute the stimulus waveforms
%        * Computation and broadcasting info about # conditions, total
%          stimulus duration, Max SPL, etc.
%
%    MakestimZWSPL renders P ready for D/A conversion by adding the following 
%    fields to P
%           Fsam: sample rate [Hz] of all waveforms. This value is
%                 determined by carrier & modulation freqs, but also by
%                 the Experiment definition P.Experiment, which may 
%                 prescribe a minimum sample rate needed for ADC.
%           Fcmp: frequencies [Hz] of all the components in a column array
%            SPL: SPLs [ms] of all the conditions
%       Waveform: Waveform object array containing the samples in SeqPlay
%                 format.
%    Attenuation: scaling factors and analog attuater settings for D/A
%   Presentation: struct containing detailed info on stimulus order,
%                 broadcasting of D/A progress, etc.
% 
%   See also Waveform/maxSPL, Waveform/play, sortConditions, 


P2 = []; % a premature return will result in []
if isempty(P), return; end
figh = P.handle.GUIfig;
EXP = P.Experiment;

% check & convert params. Note that helpers like maxSPL
% report any problems to the GUI and return [] or false in case of problems.

%     i_________ToneComplex_______: '_____ToneComplex_____'
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
ReportFsam(figh, P.Fsam/1e3); % in kHz
% ---SPL 
[P.SPL, Mess] = EvalStepper(P.StartSPL, P.StepSPL, P.EndSPL, 'linear');
if isequal('nofit', Mess),
    GUImessage(figh, {'No integer # SPL steps fit in bounds. Change tolerance or SPL step specs.'}, 'error', ...
        {'SPLstart' 'SPLstep' 'SPLend' });
elseif isequal('largstep', Mess),
    GUImessage(figh, {'Stepsize too big to fit in bounds.'}, 'error', ...
        {'SPLstart' 'SPLstep' 'SPLend' });
end
if ~isempty(Mess), return; end
Ncond = size(P.SPL,1);
if Ncond>100,
    GUImessage(figh, {'Too many (>100) SPL values. Increase step or decrease range.'}, 'error', ...
        {'SPLstart' 'SPLstep' 'SPLend' });
    return;
end
%   i_________Dur_______: '_____Dur_____'
if ~EvalDurPanel(figh, P, Ncond), 
    return;
end
% Compute the waveforms (P.Waveform) 
P = local_waveform(P, F_reg, Freq_AbsTol);

% report duration
totBaseline = sum(SameSize(P.Baseline,[1 1])); % sum of pre- & post-stim baselines
Ttotal=ReportPlayTime(figh, Ncond, P.Nrep, P.ISI, totBaseline);

% presentation order etc
P = sortConditions(P, 'SPL', 'SPL' , 'dB SPL', 'Lin');

% Levels and active channels (must be called *after* adding the baseline waveforms)
[mxSPL P.Attenuation] = maxSPL(P.Waveform, P.Experiment);
okay = CheckSPL(figh, P.SPL, mxSPL, [], '', {'StartSPL' 'EndSPL'});
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
Ncond = size(P.SPL,1);
P.StartPhase = rand(1,Ncmp); % starting phases of zwuis components
AllChan = channelSelect(P.DAC, 'LR');
chanDelay = ITD2delay(P.ITD, P.Experiment); % ITD-realizing delays for indiv D/A channels
for DAchan=AllChan,
    ichan = 1+double(P.DAC(1)=='B' && DAchan=='R');
    for icond = 1:Ncond,
        [chdel, SPL, risedur, falldur] = channelSelect(DAchan, chanDelay, P.SPL(icond,:), P.RiseDur, P.FallDur);
        ITDphaseShift = -1e-3*chdel*P.Fcmp; % delay -> phase in cycles 
        [DL, Dphi] = calibrate(P.Experiment, P.Fsam, DAchan, P.Fcmp);
        LinAmp = dB2A(SPL+DL)*sqrt(2); % calibrated linear amplitude
        ph0 = P.StartPhase(:)+Dphi(:)+ITDphaseShift; % calibrated starting phase, including ITD
        % function T = tonecomplex(Amp, Freq, Ph0, Fsam, Dur);
        w = tonecomplex(LinAmp, P.Fcmp, ph0, P.Fsam, P.BurstDur+dt); % one extra sample to avoid compulsive complaints by exactGate
        %function X = ExactGate(X, Fsam, BurstDur, Tstart, RiseDur, FallDur);
        w = ExactGate(w, P.Fsam, P.BurstDur, 0, P.RiseDur, P.FallDur);
        Param = CollectInStruct(LinAmp, ph0);
        P.Waveform(icond, ichan) = Waveform(P.Fsam, DAchan, max(abs(w)), SPL, Param, {w}, 1);
        % trailing zeros
        P.Waveform(icond, ichan) = AppendSilence(P.Waveform(icond, ichan), P.ISI);
    end
end

function  [Rise, Fall] = local_ramps(w, Nsam);
Wrise = sin(linspace(0,pi/2,Nsam).').^2;
Wfall = flipud(Wrise);
Rise = Wrise.*w(end+1+(-Nsam:-1));
Fall = Wfall.*w(end+1+(-Nsam:-1));
 





