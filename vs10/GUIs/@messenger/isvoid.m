function iv = isvoid(M);
% messenger/isvoid - true for void messenger
%    isvoid(M) returns True if M is a void messenger, false otherwise.

iv = isempty(M.Name);



