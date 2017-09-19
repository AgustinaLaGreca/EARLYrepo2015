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
    disp('Void CycleList object');
    return;
end

disp(['CycleList  ''' Q.Name ''' with callback ''' char(Q.Callback) '''.']);
for ii=1:numel(Q.Items),
    disp(['   item ' num2str(ii) ': ' Q.Items(ii).Label]);
end






