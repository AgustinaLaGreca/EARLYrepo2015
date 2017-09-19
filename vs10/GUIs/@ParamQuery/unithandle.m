function h = unithandle(Q);
% ParamQuery/unithandle - handle of unit control of rendered paramquery
%   unithandle(Q) returns the handle of the unit control of paramQuery Q. 
%   If Q has no unit control, [] is returned. Q must be rendered. 
%
%   For arrays Q, unithandle(Q) equals [unithandle(Q(1)), unithandle(Q(2)) .. ].
%   Note that this array contains less elements than Q if any of the Q(k)
%   contains no unit control.
%
%   See also ParamQuery/draw.

h = [];
for ii=1:numel(Q),
    q = Q(ii);
    if isempty(q.uiHandles),
        error(['Paramquery ''' q.Name ''' is not currently rendered.']);
    end
    h = [h getFieldOrDefault(q.uiHandles, 'Unit', [])]; 
end
    
    


