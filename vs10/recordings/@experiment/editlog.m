function editlog(E, Txt);
% experiment/editlog - edit experiment log file
%    editlog(E) opens E's log file in the Mtlab editor.
%
%    See also experiment/addtolog, experiment/insertnote,
%    experiment/status, textWrite.

if isvoid(E),
    error('Cannot type log of void experiment.');
end
edit([filename(E) '.log']);












