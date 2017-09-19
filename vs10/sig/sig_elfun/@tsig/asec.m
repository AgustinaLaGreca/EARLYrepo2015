function P = asec(P, varargin);
% tsig/asec - ASEC for tsig objects.
% 
% See also tsig.
% 
for ii=1:nchan(P),
    P.Waveform{ii} = feval(@asec, P.Waveform{ii}, varargin{:});
end;
 
