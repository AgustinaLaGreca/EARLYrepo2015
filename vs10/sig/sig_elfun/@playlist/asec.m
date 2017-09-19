function P = asec(P, varargin);
% playlist/asec - ASEC for playlist objects.
% 
% See also playlist.
% 
for ii=1:nwave(P),
    P.Waveform{ii} = feval(@asec, P.Waveform{ii}, varargin{:});
end;
 
