function P = rem(P, varargin);
% playlist/rem - REM for playlist objects.
% 
% See also playlist.
% 
for ii=1:nwave(P),
    P.Waveform{ii} = feval(@rem, P.Waveform{ii}, varargin{:});
end;
 
