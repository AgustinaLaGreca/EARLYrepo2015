function P = cot(P, varargin);
% tsig/cot - COT for tsig objects.
% 
% See also tsig.
% 
for ii=1:nchan(P),
    P.Waveform{ii} = feval(@cot, P.Waveform{ii}, varargin{:});
end;
 
