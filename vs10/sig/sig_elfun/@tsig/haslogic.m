function IR = haslogic(P);
% tsig/haslogic - true for logical channels of tsig object.
%   H=haslogic(T) is a logical row vector with H(k) is true (1) if the
%   k-th channel of T, T(k), is logical, false (0) it is numeric. 
% 
% See also tsig/isreal.
 
for ii=1:nchan(P),
    IR(ii) = islogical(P.Waveform{ii});
end;
 
