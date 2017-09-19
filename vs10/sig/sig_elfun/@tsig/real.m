function P = real(P, varargin);
% tsig/real - REAL for tsig objects.
% 
% See also tsig.
% 
for ii=1:nchan(P),
    P.Waveform{ii} = feval(@real, P.Waveform{ii}, varargin{:});
end;
 
