function P = exp(P, varargin);
% tsig/exp - EXP for tsig objects.
% 
% See also tsig.
% 
for ii=1:nchan(P),
    P.Waveform{ii} = feval(@exp, P.Waveform{ii}, varargin{:});
end;
 
