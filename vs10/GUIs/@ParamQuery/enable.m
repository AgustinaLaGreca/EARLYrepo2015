function enable(A,m)
% ParamQuery/enable - enable/disable rendered ParamQuery
%    enable(Q) or enable(Q,1) sets the enable mode of ParamQuery Q to On.
%    Enable(Q,0) sets the mode to Off.
%
%    For Arrays Q, enable(Q), enable(Q,0) and enable(Q,1) enables/disables
%    all queries. Enable(Q,M) where M is a logical array sets the enable 
%    state of element Q(k) to M(k).
%
%   See ParamQuery, ActionButton/enable.

if nargin<2, m=1; end;

NA = numel(A); 
if isequal(0,NA), % don't bother
    return;
end
if (NA>1) && isscalar(m), % multiplicate m
    m = repmat(m,1,NA);
end

if ~isequal(NA, numel(m)),
    error('Size of enable array E is incompatible with # action buttons in array A.')
end
    
mstr = {'inactive' 'on'};
ColorFactor = [0.8, 1]; 
for ii=1:NA,
    a = A(ii);
    if isfield(a.uiHandles, 'Edit') && ~isempty(a.String), % edit control
        hedt = a.uiHandles.Edit; % handle to Edit field
        CLR = getFieldOrDefault(A(ii).uicontrolProps, 'BackgroundColor', [1 1 1]); % def color is white
        set(hedt,'BackgroundColor', ColorFactor(m(ii)+1)*CLR, 'enable', mstr{m(ii)+1});
    end
    if isfield(a.uiHandles, 'Unit') && iscell(a.Unit), % toggle button
        hbtn = a.uiHandles.Unit; % handle to Btn
        CLR = getFieldOrDefault(A(ii).uicontrolProps, 'BackgroundColor', get(A(ii).Parent, 'BackgroundColor')); % def color is that of parent fig
        set(hbtn,'BackgroundColor', ColorFactor(m(ii)+1)*CLR, 'enable', mstr{m(ii)+1});
    end
    handles = struct2cell(A(ii).uiHandles);
end






