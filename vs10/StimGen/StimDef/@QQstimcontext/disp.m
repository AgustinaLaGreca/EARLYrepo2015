function disp(Q)
% stimContext/disp - DISP for ActionButton objects.
%
%   See ParamQuery.

if numel(Q)>1,
    for ii=1:numel(Q),
        disp(Q(ii));
    end
    return
end
% single element from here

if isempty(Q.DAchan),
    disp('Void stimcontext object');
    return;
end

disp('----stimContext----')
disp(struct(Q));
disp('-------------------')






