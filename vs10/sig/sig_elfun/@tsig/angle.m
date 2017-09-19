function P = angle(P, varargin);
% tsig/angle - ANGLE for tsig objects.
% 
% See also tsig.
% 
for ii=1:nchan(P),
    P.Waveform{ii} = feval(@angle, P.Waveform{ii}, varargin{:});
end;
 
