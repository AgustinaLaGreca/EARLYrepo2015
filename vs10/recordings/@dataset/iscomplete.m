function ic=iscomplete(D);
% Dataset/isvoid - true for complete (uninterrupted) dataset
%   iscomplete(D) returns TRUE if D is a complete dataset, FALSE otherwise.
%   A complete dataset is a dataset whose recording has not been
%   interrupted, resulting in the presence of recordings for all specified 
%   stimulus conditions and repetitions.
%
%   See Dataset/NrepRec, Dataset/NpresRecorded.

ic = isempty(D.Stopped);



