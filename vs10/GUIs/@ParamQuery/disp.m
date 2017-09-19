function disp(Q)
% ParamQuery/disp - DISP for ParamQuery objects.
%
%   See ParamQuery.

if numel(Q)>1,
    for ii=1:numel(Q),
        disp(Q(ii));
    end
    return
end

% single element from here
if isempty(Q),
    [M,N] = size(Q);
    disp([num2str(M) ' x ' num2str(N) ' ParamQuery object']);
    return
elseif isvoid(Q),
    disp('Void ParamQuery object');
    return;
end

if isempty(Q.Unit),
    unitstr = '';
elseif ischar(Q.Unit),
    unitstr = ['(' Q.Unit ')'] ;
else,
    unitstr = ['(' Q.Unit{1}];
    for ii=2:length(Q.Unit),
        unitstr = [unitstr '|' Q.Unit{ii}];
    end
    unitstr = [unitstr ')'];
end
disp(['ParamQuery  ' Q.Name ' ' unitstr]);







