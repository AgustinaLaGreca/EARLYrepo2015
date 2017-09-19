function P = csc(P, varargin);
% tsig/csc - CSC for tsig objects.
% 
% See also tsig.
% 
for ii=1:nchan(P),
    P.Waveform{ii} = feval(@csc, P.Waveform{ii}, varargin{:});
end;
 
