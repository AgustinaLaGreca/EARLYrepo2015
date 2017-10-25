function [EditStr, UnitStr, Enable]=grab(Q);
% ParamQuery/grab - get string values from rendered paramquery.
%   [EditStr, UnitStr, Enab]=grab(Q) gets the string properties of of the edit
%   uicontrol and, if present, the unit toggle and returns them in EditStr
%   and UnitStr, respectively. Blanks are stripped off edit strings.
%   Fixed ('read-only') unit strings and disabled toggle buttons result in
%   UnitStr==''.  Absent edit controls result in EditStr==''.
%   Enab is True when the edit field is enabled, False if not.
%
%   Unlike paramquery/read, paramquery/grab does not interpret any values 
%   nor checks their range of validity.
%
%   See also ParamQuery/read, GUIgrab.

if isempty(Q.uiHandles) || ~isfield(Q.uiHandles, 'Unit') ...
        || ~ishandle(Q.uiHandles.Unit),
    error(['Paramquery ''' Q.Name ''' is not currently rendered in a GUI.']);
end

[EditStr, UnitStr, Enable] = deal('', '', true); % defaults

% Only return non-empty EditStr if there is a edit control, i.e. if
% Q.String is non-empty (see ParamQuery and see ParamQuery/draw).
if ~isempty(Q.String),
    EditStr = get(Q.uiHandles.Edit, 'string');
    EditStr = trimspace(EditStr);
    Enable = get(Q.uiHandles.Edit, 'Enable');
end

% only return non-empty UnitStr of unit is controled by a toggle, i.e.,
% when Q.Unit is a cell array of strings.
if iscellstr(Q.Unit), 
%     enableState = get(Q.uiHandles.Unit,'enable'); % toggle object
%     if strcmp(enableState,'on') %|| strcmp(enableState,'inactive')% only return String if toggle is enabled
        UnitStr = get(Q.uiHandles.Unit, 'string');
%     end
end









