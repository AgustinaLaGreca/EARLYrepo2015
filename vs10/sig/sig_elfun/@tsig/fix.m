function P = fix(P, varargin);
% tsig/fix - FIX for tsig objects.
% 
% See also tsig.
% 
for ii=1:nchan(P),
    P.Waveform{ii} = feval(@fix, P.Waveform{ii}, varargin{:});
end;
 
