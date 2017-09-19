function P = sinh(P, varargin);
% playlist/sinh - SINH for playlist objects.
% 
% See also playlist.
% 
for ii=1:nwave(P),
    P.Waveform{ii} = feval(@sinh, P.Waveform{ii}, varargin{:});
end;
 
