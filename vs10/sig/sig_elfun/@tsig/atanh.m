function P = atanh(P, varargin);
% tsig/atanh - ATANH for tsig objects.
% 
% See also tsig.
% 
for ii=1:nchan(P),
    P.Waveform{ii} = feval(@atanh, P.Waveform{ii}, varargin{:});
end;
 
