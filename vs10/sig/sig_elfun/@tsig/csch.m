function P = csch(P, varargin);
% tsig/csch - CSCH for tsig objects.
% 
% See also tsig.
% 
for ii=1:nchan(P),
    P.Waveform{ii} = feval(@csch, P.Waveform{ii}, varargin{:});
end;
 
