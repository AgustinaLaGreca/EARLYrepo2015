function P = mod(P, varargin);
% tsig/mod - MOD for tsig objects.
% 
% See also tsig.
% 
for ii=1:nchan(P),
    P.Waveform{ii} = feval(@mod, P.Waveform{ii}, varargin{:});
end;
 
