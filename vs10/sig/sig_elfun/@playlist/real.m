function P = real(P, varargin);
% playlist/real - REAL for playlist objects.
% 
% See also playlist.
% 
for ii=1:nwave(P),
    P.Waveform{ii} = feval(@real, P.Waveform{ii}, varargin{:});
end;
 
