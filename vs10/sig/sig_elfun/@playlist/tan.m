function P = tan(P, varargin);
% playlist/tan - TAN for playlist objects.
% 
% See also playlist.
% 
for ii=1:nwave(P),
    P.Waveform{ii} = feval(@tan, P.Waveform{ii}, varargin{:});
end;
 
