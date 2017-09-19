function P = csch(P, varargin);
% playlist/csch - CSCH for playlist objects.
% 
% See also playlist.
% 
for ii=1:nwave(P),
    P.Waveform{ii} = feval(@csch, P.Waveform{ii}, varargin{:});
end;
 
