function [W, S] = width(Q);
% ParamQuery/width - estimated width in pixels of paramquery object.
%   W=width(Q) returns the estimated width W in pixels of param query Q.
%
%   [W,S]=width(Q) also returns the widths of the individual components of
%   the query in struct S with fields Prompt, Edit, Unit.
%
%   For arrays, width(Q) returns an array holding the widths of the elements.
%
%   See also ParamQuery/height.

if numel(Q)>1,
    for ii=1:numel(Q);
        W(ii) = width(Q(ii));
    end
    W = reshape(W,size(Q));
    return;
end

% ------single Q from here------
drawEdit = ~isempty(Q.String); % do we need an Edit uicontrol?

Xpr = StringExtent([Q.Prompt ' '], Q.uicontrolProps);
Xed = StringExtent([Q.String '   '], Q.uicontrolProps);
if ~drawEdit, Xed=[0 0]; end
un = Q.Unit;
if ischar(un),
    Xun = StringExtent([un ' '], Q.uicontrolProps);
elseif iscellstr(un), % width is width of largest string in cellstr + margin needed for button rims
    Xun = StringExtent(un, Q.uicontrolProps) + StringExtent('  ', Q.uicontrolProps);
else,
    error('Q.Unit is neither char nor cell string?!');
end

W = Xpr(1)+Xed(1)+Xun(1);
S = struct('Prompt', Xpr(1), 'Edit', Xed(1), 'Unit', Xun(1));

%  Name: 'StartFrequency'
%         Prompt: 'Start:'
%         String: '12000 12001'
%           Unit: 'Hz'
%     Constraint: 'rreal/positive'
%        Tooltip: 'First frequency of series.'
%         MaxNum: 2







