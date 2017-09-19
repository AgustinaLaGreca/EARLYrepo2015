function M=GUImessenger(figh,Name);
% GUImessenger - find named messenger object in GUI
%    GUImessenger(figh, 'Foo') returns the rendered messenger named Foo if
%    it exists in the GUI with. An error occurs if Foo does not exist.
%
%    See also messenger/report, GUImessage.

M = getGUIdata(figh,'Messenger');
M = M(Name);

