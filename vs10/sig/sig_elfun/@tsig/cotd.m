function P = cotd(P, varargin);
% tsig/cotd - COTD for tsig objects.
% 
% See also tsig.
% 
for ii=1:nchan(P),
    P.Waveform{ii} = feval(@cotd, P.Waveform{ii}, varargin{:});
end;
 
