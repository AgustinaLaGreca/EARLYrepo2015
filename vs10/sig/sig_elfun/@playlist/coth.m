function P = coth(P, varargin);
% playlist/coth - COTH for playlist objects.
% 
% See also playlist.
% 
for ii=1:nwave(P),
    P.Waveform{ii} = feval(@coth, P.Waveform{ii}, varargin{:});
end;
 
