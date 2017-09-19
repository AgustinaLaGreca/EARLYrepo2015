function A = setXdata(A,D)
% ActionButton/setXdata - store data asociated with an actionbutton
%    A = setXdata(A, Data) stores Data "in" the action button A. At
%    callback time, these data can be retrieved using getXdata.
%
%    A may be an array, in which case the same data D are associated with
%    each of the buttons A(k).
%
%    This is how to retrieve the data D within the callback function:
%      B = get(gcbo, 'userdata'); % retrieve action button itself
%      D = getXdata(B);
%
%   Further information on callbacks of actionbutton, right-clicking, etc,
%   is described in the help text of ActionButton.
%
%   See ActionButton, ActionButton/getXdata.

if nargout<1, error('Too few output args. Syntax is A=setXdata(A,D)'); end
if nargin<2, D=[]; end

for ii=1:numel(A),
    A(ii).Xdata = D;
end









