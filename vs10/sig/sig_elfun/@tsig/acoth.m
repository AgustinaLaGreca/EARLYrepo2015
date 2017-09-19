function P = acoth(P, varargin);
% tsig/acoth - ACOTH for tsig objects.
% 
% See also tsig.
% 
for ii=1:nchan(P),
    P.Waveform{ii} = feval(@acoth, P.Waveform{ii}, varargin{:});
end;
 
