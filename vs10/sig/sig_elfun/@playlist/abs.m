function P = abs(P, varargin);
% playlist/abs - ABS for playlist objects.
% 
% See also playlist.
% 
for ii=1:nwave(P),
    P.Waveform{ii} = feval(@abs, P.Waveform{ii}, varargin{:});
end;
 
