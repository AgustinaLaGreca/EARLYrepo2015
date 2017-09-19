function ZZ= ztrf(T, Freq, Z);
% transfer/ztrf - complex transfer function of Transfer object
%    ztrf(T) returns an array containing the complex transfer function of
%    Transfer object T.
%
%    T = ztrf(T, Freq, Z) fills T with frequency Freq (Hz) and complex
%    transfer function Z.
%
%    See Transfer, Transfer/measure, Transfer/Frequency.

if nargin==1, % get
    if ~isfilled(T),
        error('Transfer object T is not filled.');
    end
    ZZ = T.Ztrf;
else, % set
    ZZ = T;
    ZZ.Freq = Freq(:);
    ZZ.Ztrf = Z(:);
end





