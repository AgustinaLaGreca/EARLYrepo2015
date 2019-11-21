function ZP = zwuisana(St, R, DiscardBaseline, IcycleSelect, Flag);
%   zwuisana- basic zwuis analysis
%     zwuisana(St, R), with St the ZW or DZW stimulus definition as created
%     by makestimZW and R a column array of recorded values, returns a
%     struct containing the basic analysis, including discarding of ramps,
%     loopmean, fft and analysis of Rayleigh significance.
%
%     If R is a eventdata object (e.g. output of dataset/digichan), it is
%     replaced by a unary representations (zeros and ones) of the spikes
%     contained in it. 
%
%     zwuisana(St, R, DiscardBaseline) specifies whether R includes a
%     baseline portion that must be discarded or not. Default is
%     DiscardBaseline - true, as in OCT recordings.
%
%     zwuisana(St, R, DiscardBaseline, IcycleSelect) specifies which zwuis
%     cycles to use. Default is 0, which equals to using the maximum number
%     of available full zwuis cycles. A typical zwuis stimulus consists of
%     multiple reptitions of a pattern or "cycle". IcycleSelect allows the
%     analysis of one or more cycles only.
%
%     zwuisana(St, R, DiscardBaseline, IcycleSelect, 'prepare') only does
%     the preparation stage of the zwuis analysis: removal of the ramps 
%     (and of a baseline if requested); cycle selection; loop-averaging.
%     This prep-processed version of the recording can later be fed to
%     zwuisspectrum, using arbitrary frequencies to be analyzed.
%
%     Note: to know the maximum number of full zwuis cycles in the
%     dataset, first run a zwuisana on the full waveform (IcycleSelect = 0)
%     and use the NcycleZwuis from the output.
%
%     See also makestimZW, datase/apple, mscan/apple, zwuisspec.

[DiscardBaseline, IcycleSelect, Flag] ...
    = arginDefaults('DiscardBaseline, IcycleSelect, Flag', true, 0, 'complete');

[Flag, Mess] = keywordMatch(Flag, {'complete' 'prepare'});
error(Mess);

if ~isequal('ZW', St.StimType) && ~isequal('DZW', St.StimType)  && ~isequal('ZWP', St.StimType),
    error(['Zwuis analysis analysis not defined for stimulus type ' St.StimType]);
end
if ~isnumeric(IcycleSelect)
    error('IcycleSelect needs to be a numerical vector');
end
if isempty(IcycleSelect); IcycleSelect = 0; end;

dt = 1e3/St.Fsam; % ms sample period
IsSpikeData = isa(R, 'eventdata');
if IsSpikeData,
    R = local_spikes2ana(R, St);
end
[Nsam, Nz] = size(R); % # time points & # of depths stored in S
NsamPre = St.NsamRamp; % # samples before steady state
if DiscardBaseline,
    NsamPre = NsamPre + round(St.Baseline(1)/dt);
end
if size(R,1)==1, R=R.'; end
R = R(NsamPre+1:end,:); % steady state and after
Nsam = size(R,1);
NcycleZwuis = floor(Nsam/St.NsamCycle);
R = R(1:(NcycleZwuis*St.NsamCycle),:); % largest complete # cycles contained in S

if any(IcycleSelect>NcycleZwuis)
    error(['IcycleSelect requested an index larger than NcycleZwuis (%i).'...
           'Please run a full zwuisana first and use the NcycleZwuis from the output.'], NcycleZwuis);
end
if ~all(IcycleSelect == 0)
    indexes = zeros(numel(IcycleSelect)*St.NsamCycle,1);
    for iCycle = 1:numel(IcycleSelect)
        indexCycle = iCycle-1;
        targetCycle = IcycleSelect(iCycle)-1;
        indexes(1+(indexCycle*St.NsamCycle):(indexCycle*St.NsamCycle)+St.NsamCycle) =...
            1+targetCycle*St.NsamCycle:(targetCycle*St.NsamCycle)+St.NsamCycle;
    end
    R = R(indexes);
end

R = LoopMean(R,St.NsamCycle); % zwuis-cycle average of rec
if IsSpikeData, % sum rather than average across cycles
    R = R*NcycleZwuis;
end
Dur = dt*size(R,1); % ms zwuis period
Predur = dt*NsamPre; % duration [ms] of discarded pre-steady-state part
df = 1e3/Dur; % Hz spectral spacing
if isfield(St, 'SPLjump'),
    SPL = bsxfun(@plus, St.SPL, St.SPLjump(:)); % alowing for suppression stimuli
else,
    SPL = bsxfun(@plus, St.SPL, 0*St.Fzwuis(:)); % alowing for suppression stimuli
end
baseSPL = St.SPL;
% Preparation is complete. 
if isequal('prepare', Flag),
    f2i = @(f)1+round(f/df);
    i2f = @(i)(i-1)*df;
    RampDur = dt*St.NsamRamp;
    StartPhase = St.StartPhase;
    ZP = collectInStruct(dt, R, df, Predur, Dur, RampDur, NcycleZwuis, IcycleSelect, StartPhase, SPL, f2i, i2f);
    return;
end

Zsp = rfft(R);
if IsSpikeData,
    Zsp = Zsp/Zsp(1)/2; % => abs(Z.Zps) equals vector strength
    Zsp(1) = 1;
end
omega = 2*pi*Xaxis(Zsp,df) ; % angular freq in rad/s
Fprim = St.Fzwuis; % Hz
Alpha = anarayleigh(R, dt, Fprim, dt*St.NsamRamp, 10); % rayleigh significance of phase locking to stimulus
PM = pmask(Alpha<=0.001);
iprim = 1+round(Fprim/df);
MagnRaw = A2dB(abs(Zsp(iprim,:))).'; % dB re 1 nm
% estimate "threshold displacement" giving significant data
mm = sort(denan(MagnRaw+PM));
if isempty(mm), ThrUpper_dBnm = deal(nan);
else,
    ThrUpper_dBnm = mm(1); % min displ that is signif
end
ThrLower_dBnm = replaceEmpty(max(MagnRaw(isnan(PM))), nan); % max displ that is not signif
Thr_dBnm = (ThrUpper_dBnm+ThrLower_dBnm)/2;
Magn = MagnRaw + nan; % placeholder; later handed by local_apply_ref
Sens = Magn + nan; % placeholder; later handed by local_apply_ref
PhaseRaw = bsxfun(@minus, cangle(Zsp(iprim,:)).', St.StartPhase); % phase (cycle) re stim 
Cdelay = 0;  % any nonzero Cdelay will be handled by local_apply_ref
Phase = PhaseRaw + nan; % placeholder; later handed by local_apply_ref

ZP = CollectInStruct(dt, R, df, Zsp, '-apple', ...
    Thr_dBnm, ThrLower_dBnm, ThrUpper_dBnm, ...
    Dur, Cdelay, NcycleZwuis, IcycleSelect, Fprim, MagnRaw, Magn, Sens, PhaseRaw, Phase, Alpha, PM, baseSPL, SPL);

%============
function RR = local_spikes2ana(R, St)
if ~isequal(1, St.Presentation.Ncond),
    error('conversion of spike times in zwuisana is only possible for single-condition zwuis data.');
end
dt = 1e3/St.Fsam; % ms sample period
Nsam = round(1e3*St.TotDur/dt);
RR = zeros(Nsam,1);
R = eventtimes(R) - dt*St.Presentation.Onset(2); % see dataset/spiketimes
R = 1+round(R/dt);
R = R(betwixt_incl(R,1,Nsam));
RR(R) = 1;

