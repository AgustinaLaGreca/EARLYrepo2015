function P = atand(P, varargin);
% tsig/atand - ATAND for tsig objects.
% 
% See also tsig.
% 
for ii=1:nchan(P),
    P.Waveform{ii} = feval(@atand, P.Waveform{ii}, varargin{:});
end;
 
