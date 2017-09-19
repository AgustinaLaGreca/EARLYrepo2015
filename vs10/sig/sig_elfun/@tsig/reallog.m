function P = reallog(P, varargin);
% tsig/reallog - REALLOG for tsig objects.
% 
% See also tsig.
% 
for ii=1:nchan(P),
    P.Waveform{ii} = feval(@reallog, P.Waveform{ii}, varargin{:});
end;
 
