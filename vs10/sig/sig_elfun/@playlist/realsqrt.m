function P = realsqrt(P, varargin);
% playlist/realsqrt - REALSQRT for playlist objects.
% 
% See also playlist.
% 
for ii=1:nwave(P),
    P.Waveform{ii} = feval(@realsqrt, P.Waveform{ii}, varargin{:});
end;
 
