function P = sign(P, varargin);
% tsig/sign - SIGN for tsig objects.
% 
% See also tsig.
% 
for ii=1:nchan(P),
    P.Waveform{ii} = feval(@sign, P.Waveform{ii}, varargin{:});
end;
 
