function SPx = correct_offset(SPx, nsam)
% stimPresentx/correct_offset - correction for small offset
%    SPx = correct_offset(SPx, nsam) subtracts nsam from the offsets. See
%    sortConditions for why this can be useful.
%
%   See stimpresent, sortconditions.

SPx.SamOffset = SPx.SamOffset - nsam;







