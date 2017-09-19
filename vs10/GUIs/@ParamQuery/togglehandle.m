function h = togglehandle(Q);
% ParamQuery/togglehandle - handle of toggle control of rendered paramquery
%   togglehandle(Q) returns the handle of the pushbutton uicontrol of paramQuery Q. 
%   If Q has no such uicontrol, [] is returned. Q must be rendered. 
%   A paramquery is considered a toggle if its Unit property is a cellstring
%   and its Edit property is empty.
%
%   For arrays Q, togglehandle(Q) equals [togglehandle(Q(1)), togglehandle(Q(2)) .. ].
%   Note that this array contains less elements than Q if any of the Q(k)
%   contains no toggle control.
%
%   See also Param,Query/istoggle, ParamQuery/edithandle.

h = [];
for ii=1:numel(Q),
    q = Q(ii);
    if istoggle(q),
        if isempty(q.uiHandles),
            error(['Paramquery ''' q.Name ''' is not currently rendered.']);
        end
        h = [h getFieldOrDefault(q.uiHandles, 'Unit', [])];
    end
end
    
    


