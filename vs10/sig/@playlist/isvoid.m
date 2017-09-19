function IV = isvoid(T)
% playlist/isvoid - true for void playlist object
%   isvoid(P) tests if P is a void playlist object, i.e., he result of a call
%   to playlist with no input arguments.
%
%   See playlist.

IV = isempty(T.Waveform); 

