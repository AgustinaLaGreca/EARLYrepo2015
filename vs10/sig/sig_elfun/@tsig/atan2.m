function P = atan2(P, varargin);
% tsig/atan2 - ATAN2 for tsig objects.
% 
% See also tsig.
% 
for ii=1:nchan(P),
    P.Waveform{ii} = feval(@atan2, P.Waveform{ii}, varargin{:});
end;
 
