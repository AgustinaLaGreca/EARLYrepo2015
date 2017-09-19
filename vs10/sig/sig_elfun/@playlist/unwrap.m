function P = unwrap(P, varargin);
% playlist/unwrap - UNWRAP for playlist objects.
% 
% See also playlist.
% 
for ii=1:nwave(P),
    P.Waveform{ii} = feval(@unwrap, P.Waveform{ii}, varargin{:});
end;
 
