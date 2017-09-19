function P = atan2(P, varargin);
% playlist/atan2 - ATAN2 for playlist objects.
% 
% See also playlist.
% 
for ii=1:nwave(P),
    P.Waveform{ii} = feval(@atan2, P.Waveform{ii}, varargin{:});
end;
 
