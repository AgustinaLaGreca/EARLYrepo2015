function fp=fileprefix(D);
% Dataset/fileprefix - prefix of all files related to a dataset
%   fileprefix(DS) returns a char string that is the beginning of all
%   filenames having to do with datset DS.
%
%   See Dataset/save.

EE = D.ID.Experiment;
iDataset = D.ID.iDataset;
fp = [filename(EE) '_DS' zeropadint2str(iDataset,5)];



