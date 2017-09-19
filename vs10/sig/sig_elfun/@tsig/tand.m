function P = tand(P, varargin);
% tsig/tand - TAND for tsig objects.
% 
% See also tsig.
% 
for ii=1:nchan(P),
    P.Waveform{ii} = feval(@tand, P.Waveform{ii}, varargin{:});
end;
 
