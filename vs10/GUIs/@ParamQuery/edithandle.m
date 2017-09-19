function h = edithandle(Q);
% ParamQuery/edithandle - handle of edit control of rendered paramquery
%   edithandle(Q) returns the handle of the edit control of paramQuery Q. 
%   If Q has no edit control, [] is returned. Q must be rendered. 
%
%   For arrays Q, edithandle(Q) equals [edithandle(Q(1)), edithandle(Q(2)) .. ].
%   Note that this array contains less elements than Q if any of the Q(k)
%   contains no edit control.
%
%   See also ParamQuery/draw, ParamQuery/togglehandle.

h = [];
for ii=1:numel(Q),
    q = Q(ii);
    if isempty(q.uiHandles),
        error(['Paramquery ''' q.Name ''' is not currently rendered.']);
    end
    h = [h getFieldOrDefault(q.uiHandles, 'Edit', [])]; 
end
    
    


