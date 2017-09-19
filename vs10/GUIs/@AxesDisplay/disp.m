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
    disp('Void AxesDisplay object');
    return;
end

disp(['AxesDisplay  ' Q.Name ...
    ' sized ' num2str(Q.Extent(1)) ' x ' num2str(Q.Extent(2)) ' pixels']);







