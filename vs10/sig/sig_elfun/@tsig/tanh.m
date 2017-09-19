function P = tanh(P, varargin);
% tsig/tanh - TANH for tsig objects.
% 
% See also tsig.
% 
for ii=1:nchan(P),
    P.Waveform{ii} = feval(@tanh, P.Waveform{ii}, varargin{:});
end;
 
