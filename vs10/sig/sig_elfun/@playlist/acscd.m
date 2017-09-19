function P = acscd(P, varargin);
% playlist/acscd - ACSCD for playlist objects.
% 
% See also playlist.
% 
for ii=1:nwave(P),
    P.Waveform{ii} = feval(@acscd, P.Waveform{ii}, varargin{:});
end;
 
