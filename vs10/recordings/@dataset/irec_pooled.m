function I=irec_pooled(D);
% Dataset/irec_pooled - cell array of rec numbers of pooled_dataset array
%   irec_pooled(D) or D.irec_pooled returns the recording number of 
%   pooled_dataset array D a a cell array. For a non-pooled dataset array,
%   the indices are simply returned in a cell array.
%
%   Example:
%      DS = getDS('RG12432')
%      dsa = DS({21:22 21 22}) % pooled_dataset array
%      irec(DS)
%     ans = [21 21 22] % just first irecs of pooled members
%      irec_pooled(DS)
%     ans = [1x2 double]    [21]    [22] % complete cell array of indices
%
%   See also Dataset, dataset/irec_pooled.

if ~ispooled(D),
    I = num2cell(irec(D));
else,
    for ii=1:numel(D),
        M = members(D(ii));
        I{ii} = irec(M);
    end
    I = reshape(I, size(D));
end












