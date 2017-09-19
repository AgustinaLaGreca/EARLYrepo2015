function P = realpow(P, varargin);
% playlist/realpow - REALPOW for playlist objects.
% 
% See also playlist.
% 
for ii=1:nwave(P),
    P.Waveform{ii} = feval(@realpow, P.Waveform{ii}, varargin{:});
end;
 
