function P = atan(P, varargin);
% tsig/atan - ATAN for tsig objects.
% 
% See also tsig.
% 
for ii=1:nchan(P),
    P.Waveform{ii} = feval(@atan, P.Waveform{ii}, varargin{:});
end;
 
