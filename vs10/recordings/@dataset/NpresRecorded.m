function N = NpresRecorded(DS);
% Dataset/NpresRecorded - # stimulus presentations for which recordings exist
%    NpresRecorded(DS) returns the number of stimulus presentations for 
%    which dataset DS contains recordings. For a complete dataset, this 
%    number equals Npres, the total number of stimulus presentations 
%    (i.e., S.Stim.Presentation.Npres). For interrupted datasets, however,
%    recordings typically exist for a smaller number of presentations.
%
%    See also StimPresent, Dataset/NrepRec.

if iscomplete(DS), 
    N = DS.Stim.Presentation.Npres;
else,
    N = DS.Stopped.ipres-1;
end







