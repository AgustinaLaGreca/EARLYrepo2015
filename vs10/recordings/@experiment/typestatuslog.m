function typestatuslog(E);
% experiment/typestatuslog - display text from experiment statuslog file
%    typestatuslog(E) displays the text of E's statuslog file. 
%
%    See also experiment/addststustolog, experiment/status,
%    experiment/status, textWrite.

if isvoid(E),
    error('Cannot type ststuslog of void experiment.');
end
FN = [filename(E) '.statuslog'];
if exist(FN, 'file'), 
    type(FN);
end












