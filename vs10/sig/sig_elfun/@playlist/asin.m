function P = asin(P, varargin);
% playlist/asin - ASIN for playlist objects.
% 
% See also playlist.
% 
for ii=1:nwave(P),
    P.Waveform{ii} = feval(@asin, P.Waveform{ii}, varargin{:});
end;
 
