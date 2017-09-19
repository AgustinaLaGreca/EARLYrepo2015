function P = pow2(P, varargin);
% tsig/pow2 - POW2 for tsig objects.
% 
% See also tsig.
% 
for ii=1:nchan(P),
    P.Waveform{ii} = feval(@pow2, P.Waveform{ii}, varargin{:});
end;
 
