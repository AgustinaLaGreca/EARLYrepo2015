function Vrec = APtruncate(Vrec, Vthr, dt, Twin);
% APtruncate - truncate action potentials from recorded trace
%   Vrec =  APtruncate(Vrec, Vthr, dt, Twin) surgically removes the
%   action potentials from recorded trace Vrec as follows:
%     - all local maxima exceeding Vthr are considered APs
%     - for each AP, the time of the maximum, Tmax, is determined
%     - the portion of the trace from Tmax-Twin(1) to Tmax+Twin(2) is
%       replaced by the linear interpolation between its endpoints.
%
%    Inputs:
%      Vrec: recording, If matrix, each column is truncated.
%      Vthr: threshold for AP (same units as Vrec)
%        dt: sample period of recording
%      Twin: 2-element vector [-Tbefore, Tafter] as explained before. Same
%            units as dt.

if ~isvector(Vrec) & numel(Vrec)>1, % matrix: handle columns recursively
    Y = [];
    for icol = 1:size(Vrec,2),
        Y = [Y, APtruncate(Vrec(:,icol), Vthr, dt, Twin)];
    end
    Vrec = Y;
    return;
end

% book keeping
nsam0 = round(-Twin(1)/dt); % # samples in AP before peak
nsam1 = round(Twin(2)/dt); % ditto after

% make a copy of Vrec for which all values < Vthr are fixed & small
qsmall = Vrec<Vthr;
Vp = Vrec;
Vp(qsmall) = min(Vrec(:));
% now the local maxima of Vp are the APs
imax = localmax(1, Vp); 
imax = 1+imax(:); % zero-based -> one based indexing
ievent = VectorZip(imax-nsam0, imax+nsam1); % AP start indices and AP end indices alternating
ioverlap = find(diff(ievent)<=0); % overlapping APs
ievent([ioverlap; ioverlap+1]) = []; % skip the middle man
Nspike = numel(ievent)/2;
Nsam = numel(Vrec);
for ispike=1:Nspike,
    i0 = ievent(2*ispike-1); 
    i1 = ievent(2*ispike);
    % stay within trace
    i0 = max(1,i0);
    i1 = min(i1, Nsam);
    Vrec(i0:i1) = interp1([i0 i1], Vrec([i0 i1]), i0:i1);
end



