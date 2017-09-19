function P = round(P, varargin);
% playlist/round - ROUND for playlist objects.
% 
% See also playlist.
% 
for ii=1:nwave(P),
    P.Waveform{ii} = feval(@round, P.Waveform{ii}, varargin{:});
end;
 
