function P = cosh(P, varargin);
% tsig/cosh - COSH for tsig objects.
% 
% See also tsig.
% 
for ii=1:nchan(P),
    P.Waveform{ii} = feval(@cosh, P.Waveform{ii}, varargin{:});
end;
 
