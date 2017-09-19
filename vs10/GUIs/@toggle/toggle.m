function T=toggle(h,StrArray,Str0);
% toggle - constructor for toggle object
%    toggle(h,StrArray,Str0) creates toggle object for pushbutton with
%    handle h. This function is called by BETOGGLE.


if nargin<1, % void object
    [h,StrArray,Str0, Nstr, EnableState, defaultColor] = deal([]);
    T = CollectInStruct(h,StrArray,Str0, Nstr, EnableState,defaultColor);
    T = class(T,mfilename);
    return;
end

% ---arg checking
if ~isUIcontrol(h) || ~isequal('pushbutton',get(h,'style')),
    error('Input argument h must be handle to pushbutton uicontrol.');
end
if ~iscellstr(StrArray) || isempty(StrArray),
    error('StrArray arg must be non-empty cell array of strings.');
end
if isempty(Str0),
    Str0=StrArray{1};
end
% --- store info in T & create object
Nstr = length(StrArray);
EnableState = 1;
defaultColor = get(h,'BackgroundColor');
T = CollectInStruct(h,StrArray,Str0, Nstr, EnableState, defaultColor);
T = class(T,mfilename);




