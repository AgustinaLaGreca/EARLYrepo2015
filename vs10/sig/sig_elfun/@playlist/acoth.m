function P = acoth(P, varargin);
% playlist/acoth - ACOTH for playlist objects.
% 
% See also playlist.
% 
for ii=1:nwave(P),
    P.Waveform{ii} = feval(@acoth, P.Waveform{ii}, varargin{:});
end;
 
