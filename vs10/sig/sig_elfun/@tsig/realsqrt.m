function P = realsqrt(P, varargin);
% tsig/realsqrt - REALSQRT for tsig objects.
% 
% See also tsig.
% 
for ii=1:nchan(P),
    P.Waveform{ii} = feval(@realsqrt, P.Waveform{ii}, varargin{:});
end;
 
