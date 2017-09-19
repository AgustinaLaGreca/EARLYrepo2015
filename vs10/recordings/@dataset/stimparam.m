function SP = stimparam(D)
% Dataset/stimparam - struct containing stimulus parameter of dataset(s)
%   stimparam(DS) returns the struct DS.Stim. For arrays DS, a struct array
%   is returned. This is not possible using regeular syntax [DS.Stim]!
%
%   See also Dataset/subsref.

for ii=1:numel(D),
    SP(ii) = D(ii).Stim;
end

