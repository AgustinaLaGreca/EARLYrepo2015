function P = nthroot(P, varargin);
% playlist/nthroot - NTHROOT for playlist objects.
% 
% See also playlist.
% 
for ii=1:nwave(P),
    P.Waveform{ii} = feval(@nthroot, P.Waveform{ii}, varargin{:});
end;
 
