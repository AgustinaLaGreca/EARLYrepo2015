function N = gerbilName(iGerbil);
% gerbilName - experiment name from gerbil count
%   gerbilName(k) returns the experiment name for gerbil # k
%
%   See also rawdatadir.

if isnumeric(iGerbil),
    if iGerbil<=56,
        YR = '07';
    elseif iGerbil<=90,
        YR = '08';
    elseif iGerbil<=188,
        YR = '09';
    elseif iGerbil<=inf,
        YR = '10';
    end
end

N = ['RG' YR dec2base(iGerbil,10,3)];







