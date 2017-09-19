function P = sin(P, varargin);
% tsig/sin - SIN for tsig objects.
% 
% See also tsig.
% 
for ii=1:nchan(P),
    P.Waveform{ii} = feval(@sin, P.Waveform{ii}, varargin{:});
end;
 
