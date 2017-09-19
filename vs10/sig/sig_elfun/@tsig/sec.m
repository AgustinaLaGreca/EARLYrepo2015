function P = sec(P, varargin);
% tsig/sec - SEC for tsig objects.
% 
% See also tsig.
% 
for ii=1:nchan(P),
    P.Waveform{ii} = feval(@sec, P.Waveform{ii}, varargin{:});
end;
 
