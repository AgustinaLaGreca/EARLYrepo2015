function C=emptyCyclelist(dumC, hmenuitem, CLname); %#ok<INUSL>
% CycleList/emptyCyclelist - empty rendered cycleList
%   emptyCyclelist(dumC, hmenuitem, 'Foo')  empties the cyclelist named Foo
%   rendered in the uimenu which also contains the menu item
%   having handle hmenuitem. FN is abbreviated; directory & extension are
%   stripped off. dumC is a void CycleList item serving to make addFileName
%   a CycleList method.
%
%   If the name of the cycleList is omitted, the first cycle list found is
%   used.
%
%   See StimGUiFilePulldown, CycleList, CycleList/addItem.

if nargin<3, CLname = ''; end

% first retrieve cycle list C, then empty
hP = get(hmenuitem, 'parent'); % handle to uimenu containing the cycle list (see help text)
P = get(hP,'userdata'); % userdata of Pulldown object
C = getCycleList(P, CLname);
if numel(C)>1,
    error('Multiple cycle lists found in pulldown menu. Use unique name to select the required one.');
end
% clear its contents
C = rmitem(C);










