function ic=isCycleList(X)
% isCycleList - true for cycle list
%    isCycleList(X) returns true if X is a cycle list object

ic = isequal('cyclelist', lower(class(X)));