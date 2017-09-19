function [Y, T, TM, I]=getevent(X, I, Twin);
% getevent - get single event from xspike output
%   Usage: [Y, T, TM]=getevent(X, I, Twin)
%   Inputs
%        X: struct returned by xspikes
%        I: event index (may be vector). Default is 0, meaning all events
%     Twin: time window [t0 t1] in ms re event time. Default taken from X.
%   Outputs
%     Y: waveform of event taken directly from data
%     T: corresponding time vector spanning t0..t1.
%    TM: templates, that is, medians of all events in each category. If a
%        category is  not covered selection I, its template consists of NaNs.
%
%    [Y, T, TM, I]=getevent(X, -[3 4 5]) is a shorthand for 
%    getevent(X,find(ismember([X.catIdx]), [3 4 5])), that is, a selection of
%    events based on their category. The event indices are returned in
%    array I.
%
%   See also xspikes, xspikeplot, xspikeAlign.

if nargin<2, 
   I = 0;
end
if nargin<3,
    Twin = [min(X.t0snip) max(X.t0snip)];
end
if isequal(0,I),
       I = 1:sum(X.EventCount);
elseif all(I<0), % selection of event is based on their belonging to a subset of categories
    I = find(ismember([X.catIdx], abs(I)));
end
D = getABFfromSpikeStruct(X); % retrieve ABF recording
dt = X.dt;
AllSam = cat(1,D.AD(X.isweep, X.ichan).samples); % concatenate the sweeps
%xdplot(dt, 0.9*YY,'r');
Nsam = numel(AllSam);
NsamWin = 1+round(diff(Twin)/dt);
T = linspace(Twin(1), Twin(2), NsamWin)';
Y = [];
for ii=1:numel(I),
    ie=I(ii);
    t = X.tmax(ie); % event time
    isam0 = 1+round((t+Twin(1))/dt); % first sample of event
    isam1 = isam0+NsamWin-1; % last sample
    % avoid exceeding limits
    i0 = max(1,isam0);
    i1 = min(isam1, Nsam);
    y = AllSam(i0:i1);
    %isam0, isam1, i0, i1 
    if isam0<i0, y = [ones(i0-isam0,1)*y(1); y];
    elseif isam1>i1, y = [y; ones(isam1-i1,1)*y(end)];
    end
    Y = [Y, y];
end

if nargout>2, % compute templates
    TM = nan(size(Y,1), X.Ncat);
    for icat=1:X.Ncat,
        isamcat = find(X.catIdx(I)==icat);
        if ~isempty(isamcat),
            TM(:,icat) = median(Y(:,isamcat),2);
        end
    end
end





