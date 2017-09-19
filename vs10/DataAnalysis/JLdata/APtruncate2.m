function [Vrec, tAP, Tsn, Snip, APslope] = APtruncate2(Vrec, ThrSlopes, dt, Twin, makeNan);
% APtruncate - truncate action potentials from recorded trace - slope version
%   [Vrec, tAP, Tsn, Snip, APslope] =  APtruncate2(Vrec, ThrSlope, dt, Twin, makeNan) 
%   surgically removes the action potentials from recorded trace Vrec as follows:
%     - downward slopes steeper than ThrSlope mV/ms are considered APs; if
%       ThrSlope has two elements, ThrSlope(1) is the threshold for being
%       an AP and ThrSlope(2) is the threshold for being "surgically removed"
%       (see below) without actually counting as an AP.
%     - for each AP, the time of the maximum, Tmax, is determined
%     - the portion of the trace from Tmax+Twin(1) to Tmax+Twin(2) is
%       replaced by the linear interpolation between its endpoints or, if ,
%       makeNan is true, by nans.
%     - if multiple APs overlap, i.e., Tmax1+Twin(2)>Tmax2+Twin(1),
%       their union is surgically removed (instead of removing the
%       constituent APs one by one).
%
%    Inputs
%      Vrec: recording, If matrix, each column is truncated.
%    ThrSlope: AP threshold for dVrec/dt (same units as Vrec/dt)
%            if ThrSlope has two elements, ThrSlope(1) is the 
%            threshold for being an AP and ThrSlope(2) is the 
%            threshold for being "surgically removed". Anaything in
%            between is removed (just to make sure) but not accepted as AP.
%        dt: sample period of recording
%      Twin: 2-element vector [-Tbefore, Tafter] as explained before. Same
%            units as dt.
%   makeNan: logical. If true, the APs are replaced by NaNs. If false (the
%            default), they are replaced by linear extrapolation of the
%            endpoints.
%
%    Outputs
%      Vrec: truncated recording trace
%       tAP: AP peak times in same units as dt.
%       Tsn: time [ms] array corresponding to Snip (runs from Twin(1) to Twin(2))
%      Snip: snippet matrix. Columns are AP the waveforms.
%    APslope: array holding min slope in V/S of the respective APs.


[makeNan] = arginDefaults('makeNan', 0);

if ~isvector(Vrec) && numel(Vrec)>1, % matrix: handle columns recursively
    error('Cannot handle matrices.');
end

% book keeping
nsam0 = round(-Twin(1)/dt); % # samples in AP before peak
nsam1 = round(Twin(2)/dt); % ditto after
isamSnip = (-nsam0:nsam1).'; % rel indices of snippets
Tsn = dt*isamSnip; % private time axis of snippets
APthrSlope = ThrSlopes(1); % AP threshold (the strict AP criterion)
CutThr = ThrSlopes(end); % cut-out thr (the more lenient criterion of the two)

% time derivative
%dVdt = diff(Vrec)/dt;
dVdt = diff(smoothen(Vrec, 0.05, dt)/dt); % time derivative
% replace all values of dVdt that do not meet slope criterion by max(dVdt)
qsmall = dVdt>CutThr; % use CutThr for now; only afterwards remove "APs" not meeting the APthrSlope criterion
dVdt(qsmall) = max(dVdt);
% now the AP down slopes are the local minima of dVdt
idownSlope = 1+localmax(1, -dVdt); % indices of steepest down slopes
clear dVdt % save memory
%dplot(dt, Vrec); fenceplot(dt*idownSlope, ylim, 'r');

if isempty(idownSlope), % no spikes
    tAP = [];
    Snip = [];
    APslope = [];
    return;
end

% now find the maxima in Vrec that just precede the steep down slopes.
imax = [];
for ievent=1:numel(idownSlope),
    ir = idownSlope(ievent)+(-nsam1:0); % index range of "target", i.e. chunk preceding steep slope
    if ir(1)<1 || ir(end)>numel(Vrec), continue; end % whole target must be inside Vrec
    [ipeak, Vpeak] = localmax(1,Vrec(ir)); 
    if isempty(ipeak), continue; end
    [dum ilast] = max(ipeak); % get ipeak-index of last peak
    imax = [imax, ir(1)+ipeak(ilast)];
end
imax = unique(imax);

% save snippets before cutting them out
% Snip = [];
% NsamSnip = numel(isamSnip);
% for isnip=1:numel(imax),
%     ir = imax(isnip)+isamSnip;
%     if ir(1)>0 && ir(end)<=numel(Vrec),
%         Snip = [Snip Vrec(ir)];
%     end
%     APslope(1,isnip) = min(diff(Snip(:,end)))/dt;
% end
tAP = dt*(imax-1);
tAP = tAP(:);

% surgically remove the snippets, taking special care of overlapping ones.
[Vrec, Tsn, Snip] = cutFromWaveform(dt, Vrec, tAP, Twin, makeNan);
APslope(1,:) = min(diff(Snip))/dt;

% now only apply the more strict APthrSLope criterion to the candidate APs
iok = find(APslope<=APthrSlope);
[tAP, Snip, APslope] = deal(tAP(iok), Snip(:,iok), APslope(iok));


% % surgically remove the snippets, taking special care of overlapping ones.
% ievent = vectorzip(imax-nsam0, imax+nsam1); % AP start indices and AP end indices alternating
% ioverlap = find(diff(ievent)<=0); % overlapping APs
% ievent([ioverlap; ioverlap+1]) = []; % skip the middle man
% Nspike = numel(ievent)/2;
% Nsam = numel(Vrec);
% for ispike=1:Nspike,
%     i0 = ievent(2*ispike-1); 
%     i1 = ievent(2*ispike);
%     % stay within trace
%     i0 = max(1,i0);
%     i1 = min(i1, Nsam);
%     Vrec(i0:i1) = interp1([i0 i1], Vrec([i0 i1]), i0:i1);
% end
% % xplot(tAP, Vrec, 'r'); 


