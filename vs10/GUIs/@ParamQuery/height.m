function [H, S] = height(Q);
% ParamQuery/height - estimated height in pixels of paramquery object.
%   width(Q) returns the estimated height in pixels of param query Q.
%
%   [W,S]=height(Q) also returns the heights of the individual components of
%   the query in struct S with fields Prompt, Edit, Unit.
%
%   For arrays, height(Q) returns an array holding the heights of the elements.
%
%   See also ParamQuery/width.

if numel(Q)>1,
    for ii=1:numel(Q);
        H(ii) = height(Q(ii));
    end
    H = reshape(H,size(Q));
    return;
end

Xpr = StringExtent([Q.Prompt ' '], Q.uicontrolProps);
Xed = StringExtent([Q.String '  '], Q.uicontrolProps);
un = Q.Unit;
if iscell(un), % toggle button: highest string times 1.3
    Xun = StringExtent([un{:}], Q.uicontrolProps);
else, % fixed string
    Xun = StringExtent([Q.Unit ' '], Q.uicontrolProps);
end


persistent GSP; if isempty(GSP), GSP = GUIsettings('Paramquery'); end
H = GSP.HeightFactor*max([Xpr(2),Xed(2),Xun(2)]); % "double line spacing"
S = struct('Prompt', Xpr(2), 'Edit', Xed(2), 'Unit', Xun(2));







