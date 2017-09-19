function S = horzcat(varargin)
% playlist/horzcat - HORZCAT for playlist objects.
%
%   Horzcat(P,Q, ..) or [P,Q, ..] concatenates playlists P, Q, ...
%
%   Horzcat(P,X, ..) or [P,X, ..] concatenates playlist P, numerical vector
%   X, .. . The numerical parts are incorporated as single waveforms which
%   are played once.
%
%   Note: for playlists there is no distinction between Horzcat and
%   Vertcat.
%
%   See playlist.

% all args must be playlist objects or numerical vectors
iplay = cellfun(@isplaylist, varargin);
ivect = cellfun(@isnumeric, varargin) & cellfun(@isvector, varargin);

if ~all(iplay | ivect),
    error('All input args to playsig/horzcat must be playlist objects or numerical vectors.');
end

if nargin<1, % return void playlist object per convention.
    S = playlist; 
    return;
end

S = local_playlist(varargin{1});
Nw = nwave(S);
for ii=2:nargin,
    nxt = local_playlist(varargin{ii});
    S.Waveform = [S.Waveform nxt.Waveform];
    S.iplay = [S.iplay Nw+nxt.iplay];
    S.Nrep = [S.Nrep nxt.Nrep];
    Nw = nwave(S);
end

% --------------------
function P=local_playlist(P);
% convert numerical array to playlist object
if isnumeric(P),
    P = playlist(P);
    P = list(P,[1;1]); % play waveform #1 one time
end






