function [Trec, Idx] = mostRecent(T, Thist);
% mostRecent - for given instant, determine most recent instant from set
%   [Tlast, I] = mostRecent(T, Thist) returns the largest element Tlast 
%   of array Thist that is smaller than T. The index array I says
%   Tlast = Thist(I). If T coincides with an element Thist(m), then
%   this element is returned as most recent by convention. Elements of T
%   that are smaller than any elements in Thist result in Trec=NaN and I=0.
%   Empty Thist also result in Trec=NaN and I=0.
%
%   T may be an array.
%
%   See also FIND, Categorize.

Thist = Thist(:).';
if any(diff(Thist)<=0), 
    error('Elements of Thist must be strictly increasing.');
end

% treat special, trivial, cases of an empty kind.
if isempty(Thist),
    Trec = nan(size(T));
    Idx = zeros(size(T));
    return;
elseif isempty(T),
    Trec = [];
    Idx = [];
    return;
end

Idx = interp1([min(Thist(1)-1, min(T)), Thist, max(Thist(end), max(T))+1], 0:numel(Thist)+1,T);
Idx = floor(Idx);

Trec(Idx>0) = Thist(Idx(Idx>0));
Trec(Idx==0) = nan;




