function P = nthroot(P, varargin);
% tsig/nthroot - NTHROOT for tsig objects.
% 
% See also tsig.
% 
for ii=1:nchan(P),
    P.Waveform{ii} = feval(@nthroot, P.Waveform{ii}, varargin{:});
end;
 
