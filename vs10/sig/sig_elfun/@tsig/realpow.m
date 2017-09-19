function P = realpow(P, varargin);
% tsig/realpow - REALPOW for tsig objects.
% 
% See also tsig.
% 
for ii=1:nchan(P),
    P.Waveform{ii} = feval(@realpow, P.Waveform{ii}, varargin{:});
end;
 
