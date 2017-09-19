function enable(A,m)
% ActionButton/enable - enable/disable rendered actionButton
%    enable(A) or enable(A,1) sets the enable mode of action button A and
%    its associated uimenu items to On. Enable(A,0) sets the mode to Off.
%
%    For Arrays A, enable(A), enable(A,0) and enable(A,1) enables/disables
%    all buttons. Enable(A,M) where M is a logical array sets the enable 
%    state of element A(k) to M(k).
%
%
%   See ActionButton, ActionButton/highlight.

if nargin<2, m=1; end;

NA = numel(A); 
if (NA>1) && isscalar(m), % multiplicate m
    m = repmat(m,1,NA);
end

if ~isequal(NA, numel(m)),
    error('Size of enable array E is incompatible with # action buttons in array A.')
end
    
mstr = {'off' 'on'};
for ii=1:NA,
    handles = struct2cell(A(ii).uiHandles);
    set([handles{:}],'enable', mstr{m(ii)+1});
end









