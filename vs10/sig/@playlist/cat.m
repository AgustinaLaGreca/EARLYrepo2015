function S = cat(d,varargin)
% playlist/cat - CAT for playlist objects.
%
%   Cat(1,P,Q, ..) or Cat(2,P,Q, ..) concatenates playlists P,Q,.. .
%
%   Note: for playlists there is no distinction between Horzcat and
%   Vertcat.
%
%   See playlist.

S = horzcat(varargin{:});




