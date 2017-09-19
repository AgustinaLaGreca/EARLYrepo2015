function P = cot(P, varargin);
% playlist/cot - COT for playlist objects.
% 
% See also playlist.
% 
for ii=1:nwave(P),
    P.Waveform{ii} = feval(@cot, P.Waveform{ii}, varargin{:});
end;
 
