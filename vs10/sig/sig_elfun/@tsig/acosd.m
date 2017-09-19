function P = acosd(P, varargin);
% tsig/acosd - ACOSD for tsig objects.
% 
% See also tsig.
% 
for ii=1:nchan(P),
    P.Waveform{ii} = feval(@acosd, P.Waveform{ii}, varargin{:});
end;
 
