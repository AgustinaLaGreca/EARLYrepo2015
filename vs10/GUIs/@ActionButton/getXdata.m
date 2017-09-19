function D = getXdata(A,FN)
% ActionButton/setXdata - store data asociated with an actionbutton
%    getXdata(A) returns the data previously stored "in" the action button
%    A using setXdata. A may be an array of Action buttons, in which case D
%    is an array of whatever is stored in them.
%
%    getXdata(A,'Foo') only retrieves field Foo of getXdata(A) - assuming
%    it is a struct.
%
%   See ActionButton, ActionButton/setXdata.

if nargin<2, FN=''; end

if numel(A)>1,
    for ii=1:numel(A),
        D(ii) = A(ii).Xdata;
    end
else, % separate case - D may be array
    D = A.Xdata;
end

if ~isempty(FN),
    D = [D.(FN)];
end







