function P = sinh(P, varargin);
% tsig/sinh - SINH for tsig objects.
% 
% See also tsig.
% 
for ii=1:nchan(P),
    P.Waveform{ii} = feval(@sinh, P.Waveform{ii}, varargin{:});
end;
 
