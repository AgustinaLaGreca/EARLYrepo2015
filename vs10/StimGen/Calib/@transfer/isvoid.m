function iv = isvoid(T);
% transfer/isvoid - true for void Transfer object.
%    isvoid(T) returns true for a void Transfer object T, false otherwise.
%    
%    Note that a void Transfer object is not necessarily filled; see
%    Transfer/isfilled.

iv = isempty(T.Q_stim);






