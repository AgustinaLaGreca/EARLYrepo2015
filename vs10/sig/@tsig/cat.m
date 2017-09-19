function T = cat(dim, varargin)
% tsig/cat - concatenate or append tsig objects.
%   cat(1,S,T,..) or [S,T,..] is a tsig object containing the channels of
%   tsig objects S,T, .. . See tsig/horzcat for details.
%
%   cat(2,S,T,..) or [S;T;..] is a tsig object obtained by appending the
%   waveforms of S,T, in order. See tsig/vertcat for details.
%
%   See also tsig/horzcat, tsig/vertcat, tsig/subsref, tsig/repmat.

switch dim
    case 1, T = horzcat(varargin{:});
    case 2, T = vertcat(varargin{:});
    otherwise, error('Tsig objects can only be concatenated along dimensions 1 or 2.');
end

