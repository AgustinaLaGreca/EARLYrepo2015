function P = reduce_pool(P);
% pooled_dataset/reduce_pool - convert single-element pool to ordinary dataset
%    P = reduce_pool(P) returns 
%      P if P is a multi-element pool_dataset
%      P.DS if P is a single-element pool_dataset

if numel(P.DS)<2,
    P = P.DS;
end


