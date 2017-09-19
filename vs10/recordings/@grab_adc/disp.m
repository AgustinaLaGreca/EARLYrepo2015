function disp(Q)
% grabevents/disp - DISP for Grabevents objects.
%
%   See Grabevents.

try,
    Q = download(Q,'sloppy');
end

if numel(Q)>1, % recursion
    for ii=1:numel(Q),
        disp(Q(ii));
    end
    return
end

% ======single element from here=======
if isempty(Q),
    [M,N] = size(Q);
    disp([num2str(M) ' x ' num2str(N) ' Grab_adc object']);
    return
elseif isvoid(Q),
    disp('Void Grab_adc object');
    return;
end

disp(['         ---' class(Q) ' object ---'])
disp(struct(Q));






