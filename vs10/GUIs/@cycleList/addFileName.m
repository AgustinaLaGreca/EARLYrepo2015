function C=addFileName(dumC, FN, hmenuitem, CLname); %#ok<INUSL>
% CycleList/addFileName - add filename to rendered cycle list
%   addFileName(dumC, FN, hmenuitem, 'Foo')  adds filename FN to the
%   CycleList rendered in the uimenu which also contains the menu item
%   having handle hmenuitem. FN is abbreviated; directory & extension are
%   stripped off. dumC is a void CycleList item serving to make addFileName
%   a CycleList method.
%
%   If the name of the cycleList is omitted, the first cycle list found is
%   used.
%
%   See StimGUiFilePulldown, CycleList, CycleList/addItem.

if nargin<4, CLname = ''; end


if isempty(FN), return; end % apparently the user cancelled when promted for a filename
[PP FN ee] = fileparts(FN); % strip off dir & extension from file name
% first retrieve cycle list C, then add
hP = get(hmenuitem, 'parent'); % handle to uimenu containing the cycle list (see help text)
P = get(hP,'userdata'); % userdata of Pulldown object
C = getCycleList(P, CLname);
if numel(C)>1,
    error('Multiple cycle lists found in pulldown menu. Use unique name to select the required one.');
end
C = additem(C,FN);










