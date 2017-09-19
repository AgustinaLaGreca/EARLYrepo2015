function C = getCycleList(P, Name);
% PulldownMenu/getCycleList - named cyclelist of rendered PulldownMenu
%   C = getCycleList(P, 'Foo') returns cycle list named Foo of rendered 
%   PulldownMenu P. P must be rendered. C is the fully up-to-date version
%   of the CycleList as currently rendered.
%
%   C = getCycleList(P) returns all cycle lists in P.
%
%   See also PulldownMenu, PulldownMenu/additem.
if nargin<2, Name = ''; end

if ~isSingleHandle(P.Handle),
    error('getCycleList is only valid for rendered PulldownMenu objects.');
end
ihit=0;
for ii=1:numel(P.Items),
    it = P.Items{ii};
    cl =lower(class(it));
    if isequal('cyclelist',cl)
        if isempty(Name), % any cycle list is okay
            ihit = ihit + 1;
            C(ihit) = get(P.ItemHandles(ii),'userdata');
        elseif isequal(Name, it.Name),
            C = get(P.ItemHandles(ii),'userdata');
            break;
        end
    end
end
%[dum, imatch]=getitem(P,Name);








