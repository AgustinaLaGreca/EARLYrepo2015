function P = ceil(P, varargin);
% playlist/ceil - CEIL for playlist objects.
% 
% See also playlist.
% 
for ii=1:nwave(P),
    P.Waveform{ii} = feval(@ceil, P.Waveform{ii}, varargin{:});
end;
 
