function addtolog(E, Txt);
% experiment/addtolog - add text to experiment log file
%    addtolog(E, Txt) add text Txt toE's log file. For the possible formats
%    of Txt, See textwrite.
%
%    See also experiment/typelog, experiment/status, textWrite.

textwrite([filename(E) '.log'], Txt);












