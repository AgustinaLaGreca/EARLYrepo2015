function P = atan(P, varargin);
% playlist/atan - ATAN for playlist objects.
% 
% See also playlist.
% 
for ii=1:nwave(P),
    P.Waveform{ii} = feval(@atan, P.Waveform{ii}, varargin{:});
end;
 
