function P = csc(P, varargin);
% playlist/csc - CSC for playlist objects.
% 
% See also playlist.
% 
for ii=1:nwave(P),
    P.Waveform{ii} = feval(@csc, P.Waveform{ii}, varargin{:});
end;
 
