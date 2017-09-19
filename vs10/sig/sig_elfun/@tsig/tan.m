function P = tan(P, varargin);
% tsig/tan - TAN for tsig objects.
% 
% See also tsig.
% 
for ii=1:nchan(P),
    P.Waveform{ii} = feval(@tan, P.Waveform{ii}, varargin{:});
end;
 
