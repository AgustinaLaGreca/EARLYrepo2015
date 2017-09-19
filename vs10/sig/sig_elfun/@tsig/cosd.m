function P = cosd(P, varargin);
% tsig/cosd - COSD for tsig objects.
% 
% See also tsig.
% 
for ii=1:nchan(P),
    P.Waveform{ii} = feval(@cosd, P.Waveform{ii}, varargin{:});
end;
 
