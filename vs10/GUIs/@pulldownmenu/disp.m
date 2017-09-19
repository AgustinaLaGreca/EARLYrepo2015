function disp(Q)
% Messenger/disp - DISP for Messenger objects.
%
%   See ParamQuery.

if numel(Q)>1,
    for ii=1:numel(Q),
        disp(Q(ii));
    end
    return
end
% single element from here

if isempty(Q.Name),
    disp(['Void ' class(Q) ' object']);
    return;
end

disp(['PulldownMenu  ''' Q.Name ...
    ''' with label ''' Q.Label '''.']);
for ii=1:numel(Q.Items),
    it = Q.Items{ii};
    switch lower(class(it)),
        case 'struct'
            disp(['   menu item labeled ''' it.Label '''']);
        case 'pulldownmenu'
            disp(['   puldownmenu named ''' it.Name ''' and labeled''' it.Label '''']);
        case 'cyclelist'
            disp(['   Cyclelist named ''' it.Name '''']);
    end
end






