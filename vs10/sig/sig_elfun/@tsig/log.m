function P = log(P, varargin);
% tsig/log - LOG for tsig objects.
% 
% See also tsig.
% 
for ii=1:nchan(P),
    P.Waveform{ii} = feval(@log, P.Waveform{ii}, varargin{:});
end;
 
