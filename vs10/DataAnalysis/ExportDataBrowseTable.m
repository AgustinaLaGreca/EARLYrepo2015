function [Table] = ExportDataBrowseTable(Experiment)
%ExportDataBrowseTable Export the DataBrowse table out of databrowse
%   The Table from DataBrowse is exported to a struct with 4 fields
%       text: array with each row containing one dataset and its dataset
%               info
%       StimType: array of the stimulus type of each dataset
%       RecInfo: Recording info: Analog/Digital or Both
%       StimParameters: array of struct containing for each dataset a
%           struct. This struct has the following fields:
%               iRec: The recording number of the dataset =/= DatasetID
%               IDstring: The datasetID
%               Xrange: range of values for the first varied parameter of
%                   the stimulus(i.e Freauency, Intensity etc )
%               SPL: The intensity levels used in the dataset
%               Spec1: Specific parameter 1 of the stimulus
%               Spec2: Specific parameter 2 of the stimulus
%               Dur: Info on the total duration 
%               Chan: The channel(s) used: L/R/B
%
% Created By Adriaan Lambrechts 24/11/2016


[LL RR TT] = stimlist(Experiment);
Nrec = numel(TT);
RecList = [repmat(' ',[Nrec 1]) struct2char(LL) struct2char(RR)];
Table.text = RecList;
Table.StimType = TT';
Table.RecInfo = RR';
Table.StimParameters = LL';
end

