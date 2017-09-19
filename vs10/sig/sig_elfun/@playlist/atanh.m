function P = atanh(P, varargin);
% playlist/atanh - ATANH for playlist objects.
% 
% See also playlist.
% 
for ii=1:nwave(P),
    P.Waveform{ii} = feval(@atanh, P.Waveform{ii}, varargin{:});
end;
 
