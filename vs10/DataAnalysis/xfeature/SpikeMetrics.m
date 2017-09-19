function M=SpikeMetrics(X, I, Nc, Eint, APthr, maxPrePostInt, PreSpikeThr);
% SpikeMetrics - evaluate spike metrics: max EPSP rate, AP peak, etc.
%  M = SpikeMetrics(X, I, Nc, [t0 t1], , maxPrePostInt);
%             X: return struct of xspikes
%             I: selection of events (see getevent). Negative I values 
%                indicate category membership (see getevent). If I is a 
%                cell array, its elements are used for consecutive, 
%                independent analysis of the corresponding spike groups. The 
%                results are merged and sorted in the end. Such a splitting of
%                the analysis in the case of inhomogeneous groups of events
%                (typically failures vs non-failures) can improve the PCA,
%                although an increase of Nc will basically do the same.
%                The string value '>0' means all cetagories except the zero
%                one.
%            Nc: # components of PCA used for the evaluation of time
%                derivatives.
%         t0,t1: interval (ms) for evaluation of peak EPSP rate
%         APthr: AP threshold (mV). Smaller peaks will not be
%                considered to be an AP. Default value is 1.2 times the max
%                voltage during the t0..t1 interval. Specifying APthr=[]
%                also results in this default.
% maxPrePostInt: max interval (ms) from t1 for evaluation of AP peak.
%                Default maxPrePostInt= 0.5 ms.
%   PreSpikeThr: upper boundary prespike amplitude. Any values produced by 
%                the detection algorithm that exceed PreSpikeThr will
%                be set to NaN. Default is inf (=no upper boundary).
%
%  Output struct M has some bookeeping fields such as
%              Ie: indices of selected events in input struct X.
%  The metrics are in arrays corresponding to the selected events in X:
%     MaxEPSPrate: max rate of EPSP (V/s)         
%     tMaxEPSPrate: time in trace (ms) where max EPSP rate occurs
%     MinEPSPrate: (failures only) min rate of EPSPs (nan for non-failures)
%     tMinEPSPrate: time in trace (ms) where min EPSP rate occurs
%     
%     hasAP: logical array, true for events that show AP (action potential)
%            false for "failures".
%     MaxAP: peak AP amplitude (NaNs for failures)
%     MinAP: dip AP amplitude (NaNs for failures)
%     tMaxAP: time in trace (ms) of AP dip (NaNs for failures)
%
%     MinAPrate: min rate (V/s) of AP
%     tMinAPrate: time in trace (ms) where min d(AP)/dt occurs
%     
%     itvEPSP_AP: interval (ms) between EPSP max rate and peak AP (NaNs ff)
%     itvMaxMin_AP: interval between peak and dip AP (NaNs ff)
%     
%     MaxPrespike: peak amplitude (ms) of prespike
%     tMaxPrespike: time in trace (ms) where prespike peak occur
%     PreSpikePolyOrder: polynomial order used for prespike dection
%
%     See also PlotSpikeMetrics.

    
minAPFactor = 1.2; % see help text on minAP default.
PCAcrit = 0.1; % reliative deviance of 10% between waveform and PCA rep is considered a bad PCA reconstruction
PreSpikePolyOrder = 5; % polyfit order for prespike detection
MinDTevent = 0.4; % ms minimum allowed distance between events
MaxEPSPwidth = 0.4; % ms max interval between max EPSP rate & peak EPSP 

if nargin<5, APthr = []; end; % to be evaluated from the data (see help text)
if nargin<6, maxPrePostInt = 0.5; end; % ms default max pre-post interval
if nargin<7, PreSpikeThr=inf; end % default: no upper boundary for prespikes

if isequal('>0', I),
    I = -(1:X.Ncat);
elseif iscell(I), % split analysis in groups (see help text)
    for ii=1:numel(I),
        M(ii) = SpikeMetrics(X, I{ii}, Nc, Eint, APthr, ...
            maxPrePostInt, PreSpikeThr);
    end
    M = local_merge(M); 
    M = local_gethistory(M);
    return;
end

% get the events as specified by I (see help text)
[Ev, Te, dum, Ie] = getevent(X,I); Ie = Ie';
catIdx = X.catIdx(Ie).'; % event classification provided by spikes

meanEv = mean(Ev); % needed because princomp discards mean
[Coeff, Score] = princomp(Ev);
Coeff = Coeff(:,1:Nc); % restrict to requested order
REv = Score(:,1:Nc)*Coeff'; % reconstructed events
REv = REv + SameSize(meanEv, REv); % re-insert means
dt = X.dt; % sample period in ms
Etimes = X.tmax(Ie)'; % event times in row vector
Nev = numel(Ie); % total # events considered

% ----resemblance between data and PCA rep
PCAresidue = var(Ev-REv);
dataVar = var(Ev); 
poorPCA = PCAresidue >= PCAcrit*dataVar; 

% ----max EPSP rate------
iepsp = find(betwixt(Te, Eint)); % sample indices spanning EPSP range
Te_e = Te(iepsp); % idem time points re event times
REv_e = REv(iepsp, :); % idem rec events
% plot(Te, Ev(:,1:10),'linewidth',2);
% xplot(Te, REv(:,1:10));
% xplot(Te_e, REv_e(:,1:10), 'k', 'linewidth',3);
[MaxEPSPrate imax] = max(diff(REv_e)/dt); % in mV/ms = V/s
tMaxEPSPrate = dt/2+Te_e(imax)'; % time of peak EPSP rate re event time
tMaxEPSPrate = tMaxEPSPrate + Etimes; % idem, absolute time

% ----max AP--------
if isempty(APthr), % compute AP threshold from EPSP interval
    APthr = minAPFactor*max(Ev(iepsp));
end
iap = find(betwixt(Te, Eint(2)+[0 maxPrePostInt])); % sample indices spanning potential AP range
Te_ap = Te(iap); % time interval re event time for an AP to occur
Ev_ap = Ev(iap,:); % recording during that interval
[MaxAP, imax] = max(Ev_ap); % peak amplitude in that interval
tMaxAP = Te_ap(imax)' + Etimes; % time of peak
hasAP = (MaxAP>=APthr); % are the peaks APs?
MaxAP(~hasAP) = nan; % not an AP -> discard MaxAP 
tMaxAP(~hasAP) = nan; % idem tMaxAP
itvEPSP_AP = tMaxAP - tMaxEPSPrate; % interval between EPSP max rate & peak AP
% ----tInflexion, tMaxEPSP, MaxEPSP
[tInflexion, MaxEPSP, tMaxEPSP] = deal(nan(1,Nev));
for ie=1:Nev,
    if hasAP(ie),
        timespan = tMaxEPSPrate(ie) +[0 0.75]*itvEPSP_AP(ie);
        iflx = find(betwixt(Te+Etimes(ie), timespan)); % sample indices spanning inflexion time
        if numel(iflx)<2, continue; end
        REv_flx = REv(iflx, ie); % idem event
        Ev_flx = Ev(iflx, ie); % idem reconstr event
        [dum, imin] = min(diff(REv_flx));
        tInflexion(ie) = tMaxEPSPrate(ie) + (imin-0.5)*dt;
        [MaxEPSP(ie), imaxEPSP] = max(Ev_flx(1:imin));
        tMaxEPSP(ie) = timespan(1) + (imaxEPSP-1)*dt;
    else, % no AP
        timespan = tMaxEPSPrate(ie) +[0 MaxEPSPwidth];
        imep = find(betwixt(Te+Etimes(ie), timespan)); % sample indices spanning EPSP peak region
        Ev_mep = Ev(imep, ie); % idem measured event
        [MaxEPSP(ie), imaxEPSP] = max(Ev_mep);
        tMaxEPSP(ie) = timespan(1) + (imaxEPSP-1)*dt;
    end
end

% ----min AP--------
% this probes the hyperpolarization, which occurs late. So impose no later
% limit on its occurrence.
ila = find(Te>Eint(2)); % sample indices spanning everything after EPSP
Te_la = Te(ila); % idem time range
Ev_la = Ev(ila,:); % recording during that interval
[MinAP, imin] = min(Ev_la); % dip amplitude in that interval
tMinAP = Te_la(imin)' + Etimes; % time of dip
MinAP(~hasAP) = nan; % not an AP -> discard MinAP 
tMinAP(~hasAP) = nan; % idem tMaxAP
itvMinMax_AP = tMinAP - tMaxAP;

% ----MinEPSPrate and MinAPrate-----
REv_ap = REv(iap,:); % reconstructed rec during "failed AP" interval
[MinEPSPrate imin] = min(diff(REv_ap)/dt); % in mV/ms = V/s
tMinEPSPrate = dt/2+Te_ap(imin)' + Etimes;
MinAPrate = MinEPSPrate; tMinAPrate = tMinEPSPrate;
MinEPSPrate(hasAP) = nan; % what's in a name? For non-failures, we've just ...
MinAPrate(~hasAP) = nan;  % ...determined the min rate of the AP, not EPSP
tMinEPSPrate(hasAP) = nan; tMinAPrate(~hasAP) = nan;  

% ----MaxPrespike, tMaxPrespike
ipr = find(Te<Eint(1)); % sample indices spanning everything before EPSP
Te_pr = Te(ipr); % idem time range
Ev_pr = Ev(ipr,:); % recording during that interval
%Ev_pr = detrend(Ev_pr); % prespike may ride on previous event
MaxPrespike = zeros(1,Nev); tMaxPrespike = zeros(1,Nev);
for ie = 1:Nev,
    ev = Ev_pr(:,ie);
    P = polyfit(Te_pr, ev,PreSpikePolyOrder);
    rr = polyval(P,Te_pr);
    %plot(Te_pr, ev); xplot(Te_pr, rr, 'r'); pause;
    [Mrr, imax] = max(rr);
    MaxPrespike(ie) = Mrr;
    SizePrespike(ie) = Mrr-min(rr);
    tMaxPrespike(ie) = Te_pr(imax) + Etimes(ie);
end
toobig = (SizePrespike>PreSpikeThr);
SizePrespike(toobig) = nan; MaxPrespike(toobig) = nan; tMaxPrespike(toobig) = nan;
itvPrevOnset = nan(1,Nev); % time since last stim onset - to be evaluated by local_add_stimprops below
itvPrevOffset = nan(1,Nev); % time since last stim offset - to be evaluated by local_add_stimprops below
duringStim = false(1,Nev); % does EPSP occur during stimulus? - to be evaluated by local_add_stimprops below

ExpID = X.ExpID; RecID = X.RecID; icond = X.icond;
M = CollectInStruct(ExpID, RecID, icond, X, '-', ...
    I, Ie, Nc, Eint, Nev, catIdx, '-', ...
    PCAresidue, dataVar, poorPCA, PCAcrit, '-', ...
    MaxEPSPrate, tMaxEPSPrate, tInflexion, MaxEPSP, tMaxEPSP, MinEPSPrate, tMinEPSPrate, '-', ...
    hasAP, MaxAP, tMaxAP, tMinAP, MinAP, APthr, '-', ...
    MinAPrate, tMinAPrate, '-', ...
    itvEPSP_AP, itvMinMax_AP, '-', ...
    tMaxPrespike, MaxPrespike, SizePrespike, PreSpikeThr, PreSpikePolyOrder, '-', ...
    itvPrevOnset, itvPrevOffset, duringStim, '-');
M = local_gethistory(M);
M = local_add_stimprops(M,X);

%============================================
function Mm = local_merge(M); 
Mm = M(1);
Mm.I = {M.I}; Mm.Nev = sum([M.Nev]);
% merge
FN = fieldnames(M); Nfn = numel(FN);
for im=2:numel(M),
    m = M(im); % the new M to be merged
    Nev = m.Nev; % the # events m brings in
    for ifn=1:Nfn,
        fn = FN{ifn};
        newV = m.(fn);
        % check if this field is a spike metric
        isMetric(ifn) = ~isequal('Eint',fn) ...
            && (isnumeric(newV) || islogical(newV)) ...
            && isequal([1,Nev], size(newV)); 
        if isMetric(ifn),
            Mm.(fn) = [Mm.(fn) newV];
        end
    end
end
% sort according to tMaxEPSPrate
[dum, isort] = sort(Mm.tMaxEPSPrate);
for ifn=1:Nfn,
    fn = FN{ifn};
    if isMetric(ifn),
        Mm.(fn) = Mm.(fn)(isort);
    end
end

function M = local_gethistory(M); 
% for each event determine the recent history
M.itvPrevEvent = [nan diff(M.tMaxEPSPrate)]; % interval since last event
iap = [0 find(M.hasAP)];% indices of PAs preceded by a zero
APmeter = cumsum(double(M.hasAP(1:end-1))); % APmeter(k) is # APs preceeding event k+1;
M.iprevAP = [0 iap(1+APmeter)]; % index of most recent AP
itemp = M.iprevAP; itemp(itemp==0) = M.Nev+1;
tt = [M.tMaxAP nan];
M.tPrevAP = tt(itemp);
M.itvPrevAP = M.tMaxEPSPrate-M.tPrevAP; % interval between previous AP and current event
M.itvAPPrevAP = M.itvPrevAP + pmask(M.hasAP); % same metric, but restricted to APs

function M=local_add_stimprops(M,X);
% add stimulus properties
DS = sgsrdataset(M.ExpID, M.RecID);
M.qq_____STIM_____________ = '______________________';
M.DS = DS;
M.StimType = DS.StimType;
[dum M.StimInfo] = DSinfoString(DS);
M.isweep = X.isweep;
M.SweepDur = X.SweepDur;
M.OrigrepDur = DS.repdur;
M.burstDur = local_reduceifyoucan(DS.burstDur, DS.dachan, M.icond);
M.SPL = local_reduceifyoucan(getFieldOrDefault(DS, 'SPL', nan), DS.dachan, M.icond);
M.Fcar = local_reduceifyoucan(DS.carfreq, DS.dachan, M.icond);
M.Xname = DS.X.Name;
M.delay = local_reduceifyoucan(getFieldOrDefault(DS, 'delay', 0), DS.dachan, M.icond);
M.Xval = local_reduceifyoucan(DS.Xval, DS.dachan, M.icond);
M.Xunit = DS.Xunit;
if isnan(M.SPL),
    if isequal('SPL', M.StimType),
        M.SPL = M.Xval;
    else,
        warning('Cannot determine SPL of stimulus');
    end
end
[M.onset, M.offset] = deal([]);
for ii=1:numel(M.isweep),
    isw = M.isweep(ii);
    if isw>1, % stimulus starts only at 2nd sweep; first is pre-stim spont
        ons = M.delay + (isw-1)*M.SweepDur;
        ofs = ons+M.burstDur(1);
        M.onset =  + [M.onset ons];
        M.offset =  + [M.offset ofs];
        iafterons = (M.tMaxEPSPrate>ons); % events ocurring after stimulus onset
        iafterofs = (M.tMaxEPSPrate>ofs); % events ocurring after stimulus onset
        M.itvPrevOnset(iafterons) = M.tMaxEPSPrate(iafterons)-ons;
        M.itvPrevOffset(iafterofs) = M.tMaxEPSPrate(iafterofs)-ofs;
        M.duringStim(iafterons & ~iafterofs) = true;
    end
end


function X = local_reduceifyoucan(X, dachan, icond);
% reduces stumulus param to value of current condition & channel, if applicable
if dachan>0 && size(X,2)>1, % monaural stim, reduce X to that of active DA chan
    X = X(:,dachan);
end
if nargin>2 && size(X,1)>1, % X was varied; reduce to that of current condition
    X = X(icond,:);
end
% if X exists of 2 identical colms (2 channels), reduce to one
if isequal(X(:,1), X(:,end)),
    X = X(:,1);
end





