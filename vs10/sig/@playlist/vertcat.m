function S = vertcat(varargin)
% playlist/vertcat - VERTCAT for playlist objects.
%
%   Vertcat(P,Q, ..) or [P;Q; ..] concatenates playlists P, Q, ...
%
%   Note: for playlists there is no distinction between Horzcat and
%   Vertcat.
%
%   See playlist.

S = horzcat(varargin{:});




