function Yspt = SpikeTrigAv(dt, Y, Spt, Twin);
% SpikeTrigAv - spike triggered average
%   SpikeTrigAv(dt, Y, Spt, Twin) computes a spike triggered average of
%   recording Y.
%   Inputs:
%      dt: sample period 
%       Y: recording to be averaged
%     Spt: spike arrival time array in same time unit as dt
%    Twin: window [-Tpre, Tpost] relative to Spt used for averaging.
%
%    Spike times whose time windows exceed the limits of Y are ignored.
%
%    If Y is a matrix, its rows are separately spike-time-averaged.

% Y is matrix: handle columns recursively
if ~isvector(Y) && ~isempty(Y),
    Ncol = size(Y,2);
    Av = [];
    for icol=1:Ncol,
        av = SpikeTrigAv(dt, Y(:,icol), Spt, Twin);
        Av = [Av, av];
    end
    return;
end

% Y is array from here
Y = Y(:).'; % row vector 
Nsam = numel(Y);
iSpt = 1+round(Spt/dt); % indices of spikes in Y
iwin = round(Twin/dt);
iSpt = iSpt([iSpt+min(iwin)>0] & [iSpt+max(iwin)<=Nsam]); % weed out spikes whose window exceed the borders of Y
iwin = min(iwin):max(iwin); % all the indices of the window

if isempty(iSpt),
    Yspt = nan*iwin;
else,
    [iwin, iSpt] = SameSize(iwin, iSpt(:));
    Yspt = Y(iwin + iSpt); % the rows of Yspt are the "snippets" of Y
    Yspt = mean(Yspt,1);
end
if size(Y,1)>1, % Y is column vector, make Spt one too.
    Yspt = Yspt.';
end



