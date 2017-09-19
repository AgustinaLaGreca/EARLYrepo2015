function P = imag(P, varargin);
% tsig/imag - IMAG for tsig objects.
% 
% See also tsig.
% 
for ii=1:nchan(P),
    P.Waveform{ii} = feval(@imag, P.Waveform{ii}, varargin{:});
end;
 
