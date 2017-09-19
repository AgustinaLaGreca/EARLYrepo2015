function P = acos(P, varargin);
% tsig/acos - ACOS for tsig objects.
% 
% See also tsig.
% 
for ii=1:nchan(P),
    P.Waveform{ii} = feval(@acos, P.Waveform{ii}, varargin{:});
end;
 
