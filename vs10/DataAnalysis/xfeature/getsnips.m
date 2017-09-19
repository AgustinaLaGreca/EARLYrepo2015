function [SN, Tsn] = getsnips(D, dt, eT, Twin);
% getsnips - get snippets from waveform
%    [SN, Tsn] = getsnips(D, dt, eT, [Tpre, Tpost])
%    returns matrix SN whose k-th column is the snippet of waveforms cut 
%    from D starting at (eT(k)+Tpre) and ending at (eT(k)+Tpost).
%    D must be a vector holding the waveform. dt is the sample
%    period of D. eT is an array of event times (in the same time units as
%    dt, and with eT=0 corresponding to the start of the waveform D).
%    Tpre and Tpost are the offsets (re eT) of the start and
%    end times of the events to be cut, also in the same time units.
%    Tsn is a single column array holding the time axis of the the
%    snippets. Tsn runs from Tpre to Tpost, with zero corresponding to the
%    event times. 
%
%    All eT must be inside the time range of D: 0<=eT<=dt*numel(D), but it
%    is okay if part of a snippet exceeds the range. Such out-of-range
%    samples of a snippet are set to NaN.
%  
%    Example
%      [Sn, Tsn] = getsnips(D, dt, eT, [-1 1.5]) % note that Tpre<0
%      plot(Tsn, SN) % make a plot of all the snippets superimposed.
%
%    See also CutFromWaveform.

D = D(:);
eT = eT(:).';
Nsnip = numel(eT); % # requested snippets
Tpre = Twin(1); Tpost = Twin(end);
NsamSnip = round((Tpost-Tpre)/dt); % # samples per snip
D = [nan(NsamSnip,1); D; nan(NsamSnip,1)]; % provide NaN "margins" allowing out-of range cuts
Dur = dt*numel(D);
if any(eT<0) || any(eT>Dur),
    error('Event times in eT exceed the time span of the waveform D.');
end
isamstart = NsamSnip + 1+round((eT+Tpre)/dt); % sample indices of snippets' start, accounting for the nan margins introduced above
samrange = (0:NsamSnip-1).'; % zero-based index range of single snippet
isam = bsxfun(@plus, isamstart, samrange);
SN = D(isam);
Tsn = Tpre+dt*samrange;





