function N = nsam(P,k);
% playlist/nsam - number of samples of playlist items.
%
%   nsam(P) or P.nsam returns the number of samples in the items of
%   the waveforms stored in playlist P in a row vector. Note that this is
%   indpendent of the play order as specified by playlist/list.
%
%   nsam(P,K) only returns the number of samples of the Kth waveform(s).
%
%   See playlist, Playlist/list, playlist/Nwav.

N = cellfun(@length, P.Waveform);
if nargin>1,
    N = N(k);
end






