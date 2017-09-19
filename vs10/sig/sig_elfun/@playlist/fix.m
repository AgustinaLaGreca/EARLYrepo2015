function P = fix(P, varargin);
% playlist/fix - FIX for playlist objects.
% 
% See also playlist.
% 
for ii=1:nwave(P),
    P.Waveform{ii} = feval(@fix, P.Waveform{ii}, varargin{:});
end;
 
