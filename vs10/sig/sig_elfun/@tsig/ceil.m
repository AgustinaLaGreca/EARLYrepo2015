function P = ceil(P, varargin);
% tsig/ceil - CEIL for tsig objects.
% 
% See also tsig.
% 
for ii=1:nchan(P),
    P.Waveform{ii} = feval(@ceil, P.Waveform{ii}, varargin{:});
end;
 
