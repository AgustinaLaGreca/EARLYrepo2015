function P = acos(P, varargin);
% playlist/acos - ACOS for playlist objects.
% 
% See also playlist.
% 
for ii=1:nwave(P),
    P.Waveform{ii} = feval(@acos, P.Waveform{ii}, varargin{:});
end;
 
