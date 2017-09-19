function P = asinh(P, varargin);
% tsig/asinh - ASINH for tsig objects.
% 
% See also tsig.
% 
for ii=1:nchan(P),
    P.Waveform{ii} = feval(@asinh, P.Waveform{ii}, varargin{:});
end;
 
