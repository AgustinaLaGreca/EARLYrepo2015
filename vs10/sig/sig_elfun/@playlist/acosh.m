function P = acosh(P, varargin);
% playlist/acosh - ACOSH for playlist objects.
% 
% See also playlist.
% 
for ii=1:nwave(P),
    P.Waveform{ii} = feval(@acosh, P.Waveform{ii}, varargin{:});
end;
 
