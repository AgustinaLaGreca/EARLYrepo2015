function M = find(M0, figh, Name);
% messenger/find - find named messenger in GUI
%    find(messenger(), figh, 'Foo') returns the messenger named Foo
%    rendered in the GUI with handle figh. If no such messnger object
%    exists, an empty messenger is returned.

M = getGUIdata(figh, 'Messenger', pi);
if ~isequal(pi,M),
    Nam = {M.Name};
    imatch = strmatch(Name, Nam);
    if isempty(imatch),
        M = M0([]);
    else,
        M = M(imatch);
    end
end



