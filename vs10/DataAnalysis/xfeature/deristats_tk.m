function [Chunk, Raw, AP, EPSP]=deristats_tk(D, tau, ThrAPdownRate, ThrEPSP, XchunkM);
% deristats_tk - view statistics of time-derivative of recordings
%   [Chunk, Raw, AP, EPSP]=deristats_tk(D, tau, ThrAPdownRate, ThrEPSP, XchunkM)
%   INPUTS:
%          D: ABF file in TK format (see readTKABF)
%        tau: width (ms) of smoothing hann window used for time derivative.
%             Default tau=0.05 ms.
%             negative tau: don't plot...
%    ThrAPdownRate: neg rate (V/s) values below which indicate AP
%             repolarising phase of AP. Default -10 V/s.
%    ThrEPSP: threshold V for EPSP re Median. Default: 4 IQRs.
%    XchunkM: criterion for chunk extraction (pre-analysis)
%
%   OUTPUTS:
%      Chunk: elementary wave chunks
%        Raw: pointers to waveform plus elementary stats
%         AP: action potentials
%       EPSP: EPSPs

[tau, ThrAPdownRate, ThrEPSP, XchunkM] = arginDefaults('tau/ThrAPdownRate/ThrEPSP/XchunkM', []);
tau = replaceMatch(tau, [], 0.05); % ms default window width
ThrAPdownRate = replaceMatch(ThrAPdownRate, [], -10); % V/s default threshold neg rate for distinguishing APs
ThrEPSP = replaceMatch(ThrEPSP, [], []); % default: Median + 4 IQRs (realized below)
XchunkM = replaceMatch(XchunkM, [], 5); % default Xchunk factor of 5. See Xchunk

doPlot = (tau>0);
tau = abs(tau);

tau_chunk = 2*tau; % smoothing time for chunk extraction

IqrFactor = [1.5 1]*XchunkM; % chunk criterium factors passed to Xchunks.
minInflexDropRate = 1; % V/s local drop in rate to signify an inflexion

% caching
CFN = [mfilename '_' D.ExpID];
CacheParam = {D.ExpID D.RecID D.icond tau, ThrAPdownRate, ThrEPSP ...
    tau_chunk XchunkM IqrFactor minInflexDropRate , XchunkM};
Y = getcache(CFN, CacheParam);
if ~isempty(Y),
    [Chunk, Raw, AP, EPSP] = deal([], Y.Raw, Y.AP, Y.EPSP);
    [Raw.getsam, Raw.getdsam] = samgetters(D, tau);
    if doPlot, deriplot(Raw, AP, EPSP); end
    return;
end
%MaxIntPreAP = 0.8; % ms max interval preceding AP peak in which to look for event onset

[getsam, getdsam] = samgetters(D, tau); % low-budget sample getters
dt = D.dt_ms; % time spacing in ms
[sam, T] = getsam(); % whole recording & time axis in ms
dsam = getdsam(); % smoothened derivative
%sam = catsam(D,2); '=======DEBUG========'
MeanSam = mean(sam);
MedSam = median(sam); % median value of samples, serving as "baseline"
IqrSam = iqr(sam); % interquartile range, serving as "natural spread"
Nsam = numel(sam);
if isempty(ThrEPSP), ThrEPSP = 4*IqrSam; end % see help text

[ExpID, RecID, icond] = deal(D.ExpID, D.RecID, D.icond);
Raw = CollectInStruct(ExpID, RecID, icond, '-', ...
    dt, getsam, getdsam, '-', ...
    tau, tau_chunk, XchunkM, ThrAPdownRate, ThrEPSP, '-', ...
    XchunkM, IqrFactor, minInflexDropRate, '-', ...
    Nsam, MeanSam, MedSam, IqrSam);

% start chunk extraction
[Chunk, figh] = xchunks(sam, tau_chunk/D.dt_ms, IqrFactor);
Nchunk = numel(Chunk);
% add some basic fields to the chunks
isam2t = @(isam)dt*(isam-1); % quick sample index to time converter
[Chunk.t0] = DealElements(isam2t([Chunk.istart])); % start time
[Chunk.t1] = DealElements(isam2t([Chunk.iend])); % end time
[Chunk.dur] = DealElements([Chunk.t1]-[Chunk.t0]+dt); % duration
[Chunk.y0] = DealElements(sam([Chunk.istart])); % start value
[Chunk.y1] = DealElements(sam([Chunk.iend])); % end value
[Chunk.dy] = DealElements([Chunk.y1]-[Chunk.y0]); % signed span
[Chunk.isrising] = DealElements([Chunk.dy]>0);

% get more advanced metrics of each chunk.
Chunk = EvalChunks(Chunk, dsam, {'maxRate' 'tmaxRate'}, @local_maxrate, dt);
Nchunk = numel(Chunk);
% find out if a chunk has an adjacent successor; store its index.
% NOTE: these indices lose their meaning on subarrays of Chunk
hasAdjNext = [([Chunk(1:end-1).t1]==[Chunk(2:end).t0]),  false];
[Chunk.iAdjNext] = deal(nan); [Chunk(hasAdjNext).iAdjNext] = DealElements(find(hasAdjNext)+1);
% ditto predecessor
hasAdjPrev = [false, hasAdjNext(1:end-1)];
[Chunk.iAdjPrev] = deal(nan); [Chunk(hasAdjPrev).iAdjPrev] = DealElements(find(hasAdjNext)-1);


% ========APs and all that========
% determine which downward chunks are AP repolarizing phases ...
[Chunk.isAPdownFlank] = DealElements([Chunk.maxRate]<ThrAPdownRate);
% ... and which immediately precede them ...
[Chunk.isAPpreUpFlank] = DealElements([Chunk(2:end).isAPdownFlank, false] & hasAdjNext);
% ... or follow them
[Chunk.isAPpostUpFlank] = DealElements([false, Chunk(1:end-1).isAPdownFlank] & hasAdjPrev);
% inflexion point in rising pre-flank
Chunk = EvalChunks(Chunk, dsam, 'tInflex', @local_tInflex, dt, minInflexDropRate);
% Build a new type of chunks representing AP, one per chunk
iAP = find([Chunk.isAPdownFlank]); % indices of AP downward flanks in Chunk array

%%%
if ~isempty(iAP),
%%%

    AP = structpart(Chunk(iAP), {'istart' 'iend' 't0' 't1'}); % strip off other Chunk info ...
    [AP.Vpeak] = deal(Chunk(iAP).y0); % .. but re-introduce some of it while using new names.
    [AP.tpeak] = deal(Chunk(iAP).t0);
    [AP.maxDownSlope] = deal(Chunk(iAP).maxRate);
    [AP.tmaxDownSlope] = deal(Chunk(iAP).tmaxRate);
    [AP.Vdip] = deal(Chunk(iAP).y1);
    [AP.tdip] = deal(Chunk(iAP).t1);
    % add info on immediate successor of downward AP slope, if any
    hasPost = hasAdjNext(iAP); % logical indices in AP and iAP
    ipost = iAP(hasPost)+1; % indices (in Chunk) of upward post-AP slopes
    [AP(hasPost).t1] = DealElements([Chunk(ipost).t1]); % revise end time to include post-AP slope
    [AP(hasPost).iend] = DealElements([Chunk(ipost).iend]); % ditto end sample
    [AP.maxPostUpslope] = deal(nan); [AP(hasPost).maxPostUpslope] = DealElements([Chunk(ipost).maxRate]);
    [AP.tmaxPostUpslope] = deal(nan); [AP(hasPost).tmaxPostUpslope] = DealElements([Chunk(ipost).tmaxRate]);
    % add info on immediate predecessor of downward AP slope, if any
    hasPre = hasAdjPrev(iAP); % logical indices in AP and iAP
    ipre = iAP(hasPre)-1; % indices (in Chunk) of upward pre-AP slopes
    [AP(hasPre).t0] = DealElements([Chunk(ipre).t0]); % revise start time to include pre-AP slope
    [AP(hasPre).istart] = DealElements([Chunk(ipre).istart]); % ditto start sample
    [AP.tInflex] = deal(nan); [AP(hasPre).tInflex] = DealElements([Chunk(ipre).tInflex]);
    [AP.monoUpslope] = deal(true); % all APs sofar have a monotonic upward slope
    [AP.maxPreUpslope] = deal(nan); [AP(hasPre).maxPreUpslope] = DealElements([Chunk(ipre).maxRate]);
    [AP.tmaxPreUpslope] = deal(nan); [AP(hasPre).tmaxPreUpslope] = DealElements([Chunk(ipre).tmaxRate]);
    [AP.y0] = DealElements([AP.Vpeak]); [AP(hasPre).y0] = DealElements([Chunk(ipre).y0]); % start value
    % if needed, add info on short downward segment just preceding the rising AP flank
    highThr = MedSam+XchunkM*IqrSam; % high baseline
    needdownseg = (isnan([AP.tInflex]) & ([AP.y0]>highThr)); % AP-indexing logical wanting a downward segment
    tneeddownseg(1:numel(AP)) = inf; tneeddownseg(needdownseg) = [AP(needdownseg).t0]; % their t0
    hasprepre = ismember(tneeddownseg, [Chunk.t1]); % AP-indexing logical says whether downward segment is available
    [AP(hasprepre).istart] = deal(Chunk(iAP(hasprepre)-2).istart); % incorporate the downward segment by changing 1st sample idx
    [AP(hasprepre).t0] = deal(Chunk(iAP(hasprepre)-2).t0); % ditto start time
    [AP(hasprepre).tInflex] = deal(Chunk(iAP(hasprepre)-2).tmaxRate); % max downw rate is inflexion time
    [AP(hasprepre).monoUpslope] = deal(false); % these all APs do not have a monotonic upward slope
    AP = rmfield(AP,'y0');
    % again if needed, add, upward slope preceding the short downward segment
    needupseg = ~([AP.monoUpslope]); % AP-indexing logical wanting an initial upward segment
    tneedupseg(1:numel(AP)) = inf; tneedupseg(needupseg) = [AP(needupseg).t0]; % their t0
    haspreprepre = ismember(tneedupseg, [Chunk.t1]); % AP-indexing logical says whether upward segment is available
    hasprepre = hasprepre & (iAP>3); % sufficient chunks preeceeding AP must exist
    [AP(hasprepre).istart] = deal(Chunk(iAP(hasprepre)-3).istart); % incorporate the downward segment by changing 1st sample idx
    [AP(hasprepre).t0] = deal(Chunk(iAP(hasprepre)-3).t0); % ditto start time
    % get EPSPslope
    AP = EvalChunks(AP, [sam, [dsam;nan]], {'tEPSPmaxRate' 'EPSPmaxRate' 'tpeakEPSP' 'VpeakEPSP'}, @local_EPSPfromAP, dt);

%%%
else
    AP = struct('istart',NaN,'iend',NaN,'t0',NaN,'t1',NaN,'Vpeak',NaN,'tpeak',NaN,...
    'maxDownSlope',NaN,'tmaxDownSlope',NaN,'Vdip',NaN,'tdip',NaN,'maxPostUpslope',NaN,...
    'tmaxPostUpslope',NaN,'tInflex',NaN,'monoUpslope',NaN,'maxPreUpslope',NaN,...
    'tmaxPreUpslope',NaN,'tEPSPmaxRate',NaN,'EPSPmaxRate',NaN,'tpeakEPSP',NaN,...
    'VpeakEPSP',NaN,'itvAP',NaN,'itvEvent',NaN);
end;
%%%

%%%
%%%What to do with EPSP candidates, when there are now APs? Lines 158-162
%%%are problematic. But we already know there are now APs, so there can be
%%%no "swallowing".
%%%
% EPSPs & other non-AP peaks
iprepeak = find([Chunk.isrising] & ~isnan([Chunk.iAdjNext])); % CH-index array of pre-peak CHs
t0 = [Chunk(iprepeak).t0]; t1 = [Chunk(iprepeak+1).t1]; % begin & end times of peak in row vector
%%%
EPS = 1e-10;
if ~isempty(iAP),
%%%
    % find those peaks that are *completely* swallowed by APs, and remove them
    swallowed = betwixt(t0+EPS, [AP.t0].', [AP.t1].') & betwixt(t1-EPS, [AP.t0].', [AP.t1].'); % each row represents single AP
    swallowed = any(swallowed,1);  % what matters is not which AP swallows, but if any AP does it
    iup = iprepeak(~swallowed); idown = iup+1;
%%%
else
    iup = iprepeak; idown = iup+1;
end;
%%%
EPSP = struct( ...
    'istart', {Chunk(iup).istart}, 'iend', {Chunk(idown).iend}, ...
    't0', {Chunk(iup).t0}, 't1', {Chunk(idown).t1}, ...
    'dur', num2cell([Chunk(idown).t1]-[Chunk(iup).t0]), ...
    'y0', {Chunk(iup).y0}, 'y1', {Chunk(idown).y1}, ...
    'Vpeak', {Chunk(iup).y1}, ...
    'tpeak', {Chunk(iup).t1}, ...
    'Vdip', {Chunk(idown).y1}, ...
    'tdip', {Chunk(idown).t1}, ...
    'maxUpRate', {Chunk(iup).maxRate}, ...
    'tmaxUpRate', {Chunk(iup).tmaxRate}, ...
    'maxDownRate', {Chunk(idown).maxRate}, ...
    'tmaxDownRate', {Chunk(idown).tmaxRate}, ...
    'isAPpostUpFlank', {Chunk(iup).isAPpostUpFlank}, ...
    'itvAP', nan, ...
    'itvEvent', nan ...
    );
EPSP = EPSP([EPSP.Vpeak]>ThrEPSP);

%%%
%%%Here we also need a check, whether there were any APs. Because there can
%%%be itvEvents for EPSPs, but no itvAP for anything
%%%
% relative timing between events
%%%
if ~isempty(iAP),
%%%
    tAP = [AP.tpeak]; tEPSP = [EPSP.tpeak]; tEvent = sort([tAP, tEPSP]);
    [AP.itvAP] = DealElements(-EPS-mostRecent([tAP-EPS], tAP)+tAP);
    [AP.itvEvent] = DealElements(-EPS-mostRecent([tAP-EPS], tEvent)+tAP);
    if ~isempty(EPSP),
        [EPSP.itvAP] = DealElements(-EPS-mostRecent([tEPSP-EPS], tAP)+tEPSP);
    end;
%%%
else
    tEPSP = [EPSP.tpeak]; tEvent = tEPSP;
end;
%%%
if ~isempty(EPSP),
    [EPSP.itvEvent] = DealElements(-EPS-mostRecent([tEPSP-EPS], tEvent)+tEPSP);
end

% if doPlot,
%     deriplot(Raw, AP, EPSP, Chunk);
% end

% save cache, but exclude the sample getters - they take to much space
TmpGetsam = Raw.getsam;
TmpGetdsam = Raw.getdsam;
[Raw.getsam, Raw.getdsam] = deal([]);
Y = CollectInStruct(Raw, AP, EPSP);
putcache(CFN, 5e3, CacheParam, Y);
Raw.getsam = TmpGetsam;
Raw.getdsam = TmpGetdsam;

% ====== PLOT ===========


%=============================
function t = local_tInflex(dseg,c, dt, minDrop);
% evaluate time of inflection within rising flank that precedes AP falling
% flank. dseg is segment of dsam (time derivative of recording).
% MinDrop is minimum drop (V/s) in slope required for a proper inflexion.
% Nan is returned if no inflexion is detected.
t = nan;
if ~c.isAPpreUpFlank, return; end
%f6; ttt = timeaxis(dseg,dt,c.t0); plot(ttt, dseg); pause
[itmax, dsmax, isort] = localmax(1, dseg); % time in zero-based sample indices
%xplot(ttt(itmax+1), dseg(itmax+1), 'kx'); pause
if numel(itmax)<2, return; end
[itmax, dsmax] = deal(itmax(isort(1:2))+1, dsmax(isort(1:2))); % two biggest maxima; time as one-based indices
[itmax, dsmax] = sortAccord(itmax, dsmax, itmax); % chronological order
%itmax, dsmax, dsize(dseg, ttt)
%xplot(ttt(itmax), dseg(itmax), 'r*');
dseg = dseg(itmax(1):itmax(end)); % restrict segment to stuff between peaks of dseg.
%ttt = ttt(itmax(1):itmax(end)); xplot(ttt, dseg, 'r'); pause
[itmin, dsmin, isort] = localmax(1, -dseg); % minim of dseg between peaks
if isempty(itmin), return; end % no local minium detected
[itmin, dsmin] = deal(itmin(isort(1))+1, -dsmin(isort(1))); % select deepest minimum found
%xplot(ttt(itmin), dsmin, 'sk'); pause
itmin = itmin + itmax(1); % location of min re start of original, complete, dseg, as 1-based index
Drop = dsmax-dsmin; % depth of minimum re both surrounding maxima
if min(Drop)<abs(minDrop), return; end % not enough drop to call it an inflexion
t = c.t0+dt*(itmin-2); % absolute time of inflexion
%xplot(t, dsmin, 'sk', 'markerfacecolor', [1 0 0]); pause; delete(gcf);

function [maxRate, tmaxrate] = local_maxrate(dsam,c, dt);
% chunk function for extracting max abs rate. Segment cut from time derivative dsam
[maxRate, imax] = max(abs(dsam));
maxRate = sign(c.dy)*maxRate;
tmaxrate = c.t0+(imax-1)*dt;

function [t, mR, tE, mE] = local_EPSPfromAP(sam_dsam, c, dt);
% get steepest slope & its time of initial segment of AP, i.e. the EPSP part
% also get peak of EPSP and its time
sam = sam_dsam(:,1); dsam = sam_dsam(1:end-1,2);
if isnan(c.tInflex),
    t = nan;
    mR = nan;
    tE = nan;
    mE = nan;
else,
    T = dt*(c.istart:c.iend)-dt; % proper time of chunk c
    iep = find(T<c.tInflex); % sample indices of the pre-tInflex portion
    [mR, imax] = max(dsam(iep)); % steepest point of the pre-tInflex portion
    t = T(iep(imax)); % its time
    [mE, imax] = max(sam(iep)); % peak of pre-tInflex portion
    tE = T(iep(imax)); % its time
end






