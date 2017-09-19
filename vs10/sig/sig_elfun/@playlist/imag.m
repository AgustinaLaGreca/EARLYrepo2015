function P = imag(P, varargin);
% playlist/imag - IMAG for playlist objects.
% 
% See also playlist.
% 
for ii=1:nwave(P),
    P.Waveform{ii} = feval(@imag, P.Waveform{ii}, varargin{:});
end;
 
