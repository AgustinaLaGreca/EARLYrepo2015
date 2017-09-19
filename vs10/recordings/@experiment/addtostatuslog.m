function addtostatuslog(E, Txt);
% experiment/addtostatuslog - add text to experiment statuslog file
%    addtostatuslog(E, Txt) adds text Txt to E's statuslog file. For the 
%    possible formats of Txt, See textwrite.
%
%    See also experiment/typestatuslog, experiment/status, textWrite.

textwrite([filename(E) '.statuslog'], Txt);












