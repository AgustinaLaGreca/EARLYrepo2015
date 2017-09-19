function P = cos(P, varargin);
% tsig/cos - COS for tsig objects.
% 
% See also tsig.
% 
for ii=1:nchan(P),
    P.Waveform{ii} = feval(@cos, P.Waveform{ii}, varargin{:});
end;
 
