function IR = isreal(P);
% tsig/isreal - ISREAL for tsig objects.
%   isreal(T) is a logical row vector whose elements correspond to the
%   channels of T and values are true (1) for real channels, false (0) for
%   complex channels. Note: logical signals are considered real.
% 
% See also tsig/haslogic.
% 
for ii=1:nchan(P),
    IR(ii) = isreal(P.Waveform{ii});
end;
 
