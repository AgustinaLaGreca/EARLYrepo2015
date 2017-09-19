function P = sind(P, varargin);
% tsig/sind - SIND for tsig objects.
% 
% See also tsig.
% 
for ii=1:nchan(P),
    P.Waveform{ii} = feval(@sind, P.Waveform{ii}, varargin{:});
end;
 
