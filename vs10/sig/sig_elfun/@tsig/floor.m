function P = floor(P, varargin);
% tsig/floor - FLOOR for tsig objects.
% 
% See also tsig.
% 
for ii=1:nchan(P),
    P.Waveform{ii} = feval(@floor, P.Waveform{ii}, varargin{:});
end;
 
