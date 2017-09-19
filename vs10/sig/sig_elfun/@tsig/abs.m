function P = abs(P, varargin);
% tsig/abs - ABS for tsig objects.
% 
% See also tsig.
% 
for ii=1:nchan(P),
    P.Waveform{ii} = feval(@abs, P.Waveform{ii}, varargin{:});
end;
 
