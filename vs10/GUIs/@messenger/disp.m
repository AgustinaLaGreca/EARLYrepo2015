function disp(Q)
% Messenger/disp - DISP for Messenger objects.
%
%   See ParamQuery.

if numel(Q)>1,
    for ii=1:numel(Q),
        disp(Q(ii));
    end
    return
elseif numel(Q)==0,
    disp([sizeString(size(Q)) ' Messenger object'])
    return;
end
% single element from here

if isempty(Q.Name),
    disp('Void Messenger object');
    return;
end


disp(['Messenger  ' Q.Name ]);







