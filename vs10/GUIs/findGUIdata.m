function H = findGUIdata(FN, value, IncludeHidden);
% findGUIdata - find figures comtaining certain getGUIdata.
%   findGUIdata('Prop', Val) returns an array of handles of those figures
%   whose GUIdata equals Val. Prop may be "staggered" as in Prop1.field.
%
%   findGUIdata('Prop', Val, IncludeHidden) includes hidden handles if
%   IncludeHidden is true. The default is not to include hidden handles in
%   the search.
%
%   See also getGUIdata, setGUIdata.


IncludeHidden = arginDefaults('IncludeHidden', 0);
hiha = get(0, 'ShowHiddenHandles');
if IncludeHidden,
    set(0,'ShowHiddenHandles', 'on');
end
fh = findobj('type', 'figure');
set(0, 'ShowHiddenHandles', hiha);
H = [];
for ii=1:numel(fh),
    if isequal(value, getGUIdata(fh(ii), FN, now)),
        H = [H fh(ii)];
    end
end

