function P = acsc(P, varargin);
% playlist/acsc - ACSC for playlist objects.
% 
% See also playlist.
% 
for ii=1:nwave(P),
    P.Waveform{ii} = feval(@acsc, P.Waveform{ii}, varargin{:});
end;
 
