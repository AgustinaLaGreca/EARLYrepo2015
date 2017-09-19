function clearedit(Q);
% ParamQuery/clearedit - clear edit control of rendered paramquery
%   clearedit(Q) clears the edit field of the elements in Q.
%   Q must be rendered.
%
%   See also ParamQuery/read, GUIgrab.


for ii=1:numel(Q),
    q = Q(ii);
    if isempty(q.uiHandles) || ~isfield(q.uiHandles, 'Unit') ...
            || ~isSingleHandle(q.uiHandles.Unit),
        error(['Paramquery ''' q.Name ''' is not currently rendered in a GUI.']);
    end

    if isfield(q.uiHandles, 'Edit')
        set(q.uiHandles.Edit, 'string','');
    end
end

