function P = tand(P, varargin);
% playlist/tand - TAND for playlist objects.
% 
% See also playlist.
% 
for ii=1:nwave(P),
    P.Waveform{ii} = feval(@tand, P.Waveform{ii}, varargin{:});
end;
 
