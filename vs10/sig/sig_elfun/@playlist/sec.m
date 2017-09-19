function P = sec(P, varargin);
% playlist/sec - SEC for playlist objects.
% 
% See also playlist.
% 
for ii=1:nwave(P),
    P.Waveform{ii} = feval(@sec, P.Waveform{ii}, varargin{:});
end;
 
