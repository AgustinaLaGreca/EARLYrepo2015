function P = acosd(P, varargin);
% playlist/acosd - ACOSD for playlist objects.
% 
% See also playlist.
% 
for ii=1:nwave(P),
    P.Waveform{ii} = feval(@acosd, P.Waveform{ii}, varargin{:});
end;
 
