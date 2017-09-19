function B=baseline(SPx, Fld)
% stimPresentx/baseline - baseline-related values
%   B = baseline(SPx) returns struct B with fields
%         PreDur: duration (ms) of pre-stimulus baseline
%        PostDur: duration (ms) of post-stimulus baseline
%        NsamPre: # samples of pre-stimulus baseline
%       NsamPost: # samples of post-stimulus baseline
%
%   B = baseline(SPx, 'foo') only returns filed foo of this struct.
%
%   See also stimpresentx, sortConditions.

B = structpart(SPx, {'PreDur', 'PostDur',  'NsamPre', 'NsamPost'});
if nargin>2,
    B = B.(Fld);
end


