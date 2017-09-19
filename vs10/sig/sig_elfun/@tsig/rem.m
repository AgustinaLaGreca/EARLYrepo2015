function P = rem(P, varargin);
% tsig/rem - REM for tsig objects.
% 
% See also tsig.
% 
for ii=1:nchan(P),
    P.Waveform{ii} = feval(@rem, P.Waveform{ii}, varargin{:});
end;
 
