function P = acosh(P, varargin);
% tsig/acosh - ACOSH for tsig objects.
% 
% See also tsig.
% 
for ii=1:nchan(P),
    P.Waveform{ii} = feval(@acosh, P.Waveform{ii}, varargin{:});
end;
 
