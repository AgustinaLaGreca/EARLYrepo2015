function P = cos(P, varargin);
% playlist/cos - COS for playlist objects.
% 
% See also playlist.
% 
for ii=1:nwave(P),
    P.Waveform{ii} = feval(@cos, P.Waveform{ii}, varargin{:});
end;
 
