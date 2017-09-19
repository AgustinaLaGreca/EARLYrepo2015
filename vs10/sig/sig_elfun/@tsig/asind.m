function P = asind(P, varargin);
% tsig/asind - ASIND for tsig objects.
% 
% See also tsig.
% 
for ii=1:nchan(P),
    P.Waveform{ii} = feval(@asind, P.Waveform{ii}, varargin{:});
end;
 
