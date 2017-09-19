function P = cosd(P, varargin);
% playlist/cosd - COSD for playlist objects.
% 
% See also playlist.
% 
for ii=1:nwave(P),
    P.Waveform{ii} = feval(@cosd, P.Waveform{ii}, varargin{:});
end;
 
