function [S, SP] = apple(DS, Chan, iCond, ME);
% dataset/apple - analysis of ZW data
%    apple(DS, Chan, iCond) computes and displays the linear transfer from
%    analog responses to ZW (zwuis) stimuli.
%    Chan is the channel spec of the analog channel. See dataset/anachan.
%    iCond is the stimulus condition. Default iCond is 0, meaning all
%    conditions.
%
%    DS may be in array of datasets, in which case all conditions of all
%    datasets are simply "concatenated". To select different conditions
%    from different datasets, use a cell array for iCond as in
%      apple(DS(11:14), 2, {[1] [2 4 7] [2] [5]})
%    Here each cell of the last input arg selects the stimulus conditions
%    of the corresponding dataset array.
%
%    apple(DS, Chan, iCond, ME) normalizes the gain and phases using the
%    transfer function from transfer object ME.
%  
%    apple(dataset(), S, 'plot') plots a previous output S of apple and
%    returns the handles to the line objects of this plot.
%    apple(dataset(), S, 'plot', ha) uses axes ha(1) & ha(2) for this plot.
%
%    Caching: file name apple_RGxxxx.
%  
%    See also dataset/powerspec.

[iCond, ME] = arginDefaults('iCond/ME', 0, []); % default: all conditions, no normalization

if nargin>=3 && ischar(iCond) && isequal('plot', lower(iCond)) && isvoid(DS),
    S = local_plot(Chan, ME);
    return;
end

iCond = cellify(iCond); % each cell ~ element of DS
iCond = SameSize(iCond, DS); % if needed, replicate condition to match # DS elements

S = []; SP = [];
for ids=1:numel(DS),
    ds = DS(ids);
    ds = local_checkpool(ds);
    icond = iCond{ids};
    if isequal(0,icond),
        icond = 1:Ncond(ds); % icond==0 means "all conditions"
    end
    for ii=1:numel(icond),
        [s, sp] = local_apple(ds, Chan, icond(ii), ME);
        S = [S, s];
        SP = [SP, sp];
    end
end

[S, SP] = sortAccord(S, SP, [S.Fprofile] + 0.1*[S.baseSPL]+0.001*[S.SPLjump]);
if nargout<1, % plot
    local_plot(S,[]);
end

% ======================================================================
function [S, Spec] = local_apple(DS, Chan, iCond, ME); % single condition of singe DS
CP = {uniqueID(DS), datenum(time_created(DS)), Chan, iCond};
if ~isempty(ME), CP = [CP getWBdelay(ME), description(ME)]; end
CFN = ['apple_' expname(DS)];
SSp = getcache(CFN, CP);
if ~isempty(SSp), 
    [S, Spec] = deal(SSp.S, SSp.Spec);
    return; 
end
[Cspec, df, Fprim, Alpha] = local_Cspec(DS, Chan, iCond);
freq = 1e3*Xaxis(Cspec, df); % Hz
absMagn = A2dB(abs(Cspec));
absPhase = cangle(Cspec);
% from here, only bookkeeping data are extracted from DS, so reduce any
% pooled datasets.
irec_pool = irec_pooled(DS); % complete rec indices
if isa(DS, 'pooled_dataset'),
    DS = members(DS);
    DS = DS(1);
end
% undo DS's own hardware timelag
absPhase = delayPhase(absPhase, freq/1e3, timelag(anachan(DS,Chan)));
if ~isempty(ME), % apply normalization
    ctrf = eval(ME, freq,1); % complex transfer fnc; last 1 = allow out-of-range freqs
    absMagn = absMagn - A2dB(abs(ctrf));
    absPhase = absPhase - cangle(ctrf); % ctrf does not include the wideband delay of ME and it shouldn't!
end
%Freq = Xaxis(Spec,df);
isam = @(fr)1+round(1e-3*fr/df); % Hz -> 1-base sample index
iPrim = isam(Fprim); % indices of zwuis freqs
SPLprim = DS.Stim.SPL + DS.Stim.SPLjump(iCond,:); % SPL of individual zwuis components in Nzwuis row column
Gain = absMagn(iPrim,:).' - SPLprim; % Gain(icond, iprim)
Phase = bsxfun(@minus, absPhase(iPrim,:).', DS.Stim.StartPhase); % phase re stim phase
% collect all params & return struct S
ExpName = name(DS.ID.Experiment);
irec = DS.ID.iDataset;
Fprofile = DS.Stim.ProfileFreq(iCond);
ProfileType = local_proftype(DS.Stim.ProfileType);
baseSPL = DS.Stim.SPL;
SPLjump = DS.Stim.ProfileSPLjump;
if isequal([], ME),
    Reference = [];
else
    Reference = strip(ME); % strip content to save space. Call unstrip(Reference) to recover.
end
LegStr = sprintf('%2d (+%2d) dB  %0.1f kHz  %s', round(baseSPL), SPLjump, Fprofile/1e3, ProfileType);
uID = uniqueID(DS);
S = CollectInStruct(uID, ExpName, irec, irec_pool, '-zwuis', Fprim, Gain, Phase, Alpha, ...
    '-zwuisparams', Fprofile, ProfileType, baseSPL, SPLjump, LegStr, Reference);
Spec = CollectInStruct(df, Cspec);
SSp = CollectInStruct(df, S, Spec);
putcache(CFN, 200, CP, SSp);

function [Cspec, df, Fprim, RayleighAlpha, NsamSpec, NsamRamp, dt] = local_Cspec(DS, Chan, iCond);
if isa(DS, 'pooled_dataset'),
    [Cspec, df, Fprim, RayleighAlpha] = local_pooledCspec(DS, Chan, iCond);
    return;
end
bdur = DS.Stim.BurstDur(1);
[Scal, Yunit] = conversionfactor(anachan(DS, Chan));
[D, dt]= anamean(DS,Chan, iCond, 1, [0 bdur]);
NsamCyc = DS.Stim.NsamCycle; % Number of samples per zwuis base cycle
NsamRamp = DS.Stim.NsamRamp;
% discard ramps
D = Scal*D(NsamRamp+1:end-NsamRamp, :);
% round to integer # zwuis cycles & apply loopmean
NsamSpec = NsamCyc*floor(size(D,1)/NsamCyc);
D = D(1:NsamSpec,:);
Fprim = DS.Stim.Fzwuis;
% Compute mean of zwuis repetitions
D = D - mean(D(:));
D = simplegate(D,NsamRamp); % apply temporal window matching that of stimulus
[D M] = LoopMean(D,NsamCyc);
% Compute the FFT
Cspec = fft(D)*2/numel(D);
RayleighAlpha = anarayleigh(D, dt, Fprim, dt*NsamRamp, 10); % rayleigh significance of phase locking to stimulus
df = 1/(dt*size(D,1)); % frequency spacing in kHz

function [Cspec, df, Fprim, RayleighAlpha] = local_pooledCspec(P, Chan, iCond);
DS = members(P);
NP = numel(DS); 
Cspec = 0; 
for ii=1:NP,
    [csp, df, Fprim, ra, w(ii), NsamRamp, dt] = local_Cspec(DS(ii), Chan, iCond);
    Cspec = Cspec + w(ii)*csp;
end
Cspec = Cspec/sum(w);
Drec = ifft(Cspec); % reconstruct weighted sum of time signal
RayleighAlpha = anarayleigh(Drec, dt, Fprim, dt*NsamRamp, 10);

function Str = local_proftype(ProfileType);
switch ProfileType, 
    case 'single cmp', Str = '1';
    case 'low end', Str = '>';
    case 'high end', Str = '<';
    case '-', Str = '-';
end


function LH = local_plot(S, ha);
if isempty(ha),
    set(gcf,'units', 'normalized', 'position', [0.28 0.293 0.647 0.606]);
    ha(1) = subplot(2,1,1);
    ha(2) = subplot(2,1,2);
end
for ii=numel(S):-1:1,
    s = S(ii);
    PLM = pmask(s.Alpha<=0.01);
    axes(ha(1));
    freq = s.Fprim/1e3;
    LH.Gain_dot(ii) = xplot(freq, s.Gain, '.:', lico(ii), 'markersize', 4);
    LH.Gain_lin(ii) = xplot(freq, s.Gain+PLM, ['-' ploma(ii)], lico(ii), 'markersize', 5);
    xlog125([min(freq)/2 min(50, 1.5*max(freq))]);
    %legend(S.LegStr, 'location', 'west');
    axes(ha(2));
    LH.Phase_dot(ii) = xplot(freq, cunwrap(s.Phase), '.:' , lico(ii), 'markersize', 4);
    LH.Phase_lin(ii) = xplot(freq, cunwrap(s.Phase)+PLM, ['-' ploma(ii)], lico(ii), 'markersize', 5);
    grid on;
end
axes(ha(1));
%title([S(1).ExpName '  rec ' sprintf('%d ', [S.irec])]);
title([S(1).ExpName '  rec ' numcell2str([S.irec_pool])], 'interpreter', 'none');

ylabel('Gain (dB)');
legend(LH.Gain_lin, {S.LegStr}, 'location', 'northwest')
axes(ha(2));
xlabel('Frequency (kHz)');
ylabel('Phase (cycle)');
phasedelayslider(ha(2), 'kHz');
linkaxes(ha, 'x');

function ds = local_checkpool(ds);
ds = reduce_pool(ds);
ST = upper(stimtype(ds));
if ~isa(ds, 'pooled_dataset'), 
    if ~isequal('ZW', ST),
        error(['Dataset ' IDstring(ds, 'full') 'is not a ZW dataset.']);
    end
else, % pooled dataset
    if iscell(ST),
        ST = unique(ST);
        if numel(ST)>1,
            error('Not all members of pooled dataset are ZW datasets.');
        end
    end
    stref = [];
    M = members(ds);
    for ii=1:numel(M),
        st = M(ii).Stim;
        st = structpart(st, {'BaseFreq' 'ZWseed' 'Ncomp' 'Fzwuis' 'StartPhase' 'SPL' 'NsamCycle'});
        if isempty(stref),
            stref = st;
        else,
            if ~isequal(st, stref),
                st
                stref
                error('Attempting to pool incompatible datasets.');
            end
        end
    end
end


