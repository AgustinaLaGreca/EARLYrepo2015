function P = round(P, varargin);
% tsig/round - ROUND for tsig objects.
% 
% See also tsig.
% 
for ii=1:nchan(P),
    P.Waveform{ii} = feval(@round, P.Waveform{ii}, varargin{:});
end;
 
