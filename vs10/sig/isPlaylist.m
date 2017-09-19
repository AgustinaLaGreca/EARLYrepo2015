function IP = isPlaylist(P);
% isPlaylist - true for playlist object
%    isPlaylist(P) returns logical true (1) if P is a playlist object and
%    logical false (0) otherwise.
%
%    See also playlist.

IP = isa(P, 'playlist');