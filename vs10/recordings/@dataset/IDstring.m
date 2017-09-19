function IDstr = IDstring(D, flag)
% Dataset/IDstring - char string identifying a dataset within an experiment.
%   IDstring(DS) returns a string like '3-5-FS' meaning that dataset DS was
%   recorded using an FS stimulus, and that it was recording #5 obtained
%   from unit ("cell") #3.
%
%   IDstring(DS, 'full') returns more info, e.g.,
%     RG11251, Rec 4  <1-2-FS>
%
%   See also Dataset, Dataset/paramview, Dataset/stimlist.

flag = arginDefaults('flag', '');
if ~isempty(flag),
    [flag, Mess] = keywordMatch(flag, {'full','status'}, 'flag argument');
    error(Mess);
end

icStr = num2str(D.ID.iCell);
iocStr = num2str(D.ID.iRecOfCell);
idatStr = num2str(D.ID.iDataset);
if any(isnan([D.ID.iCell D.ID.iRecOfCell])) % for datagraze => dataset ID is unknown
    stat = status(D.ID.Experiment);
    icStr = num2str(stat.iCell);
    iocStr = num2str(stat.iRecOfCell+1);
    idatStr = num2str(stat.Ndataset);
end

if isfield(D.ID,'SubStimType')
    if ~(isempty(D.ID.SubStimType)|| strcmp(D.ID.SubStimType,''))
        IDstr = [icStr '-' iocStr '-' upper(stimtype(D)) '-' D.ID.SubStimType];
    else
        IDstr = [icStr '-' iocStr '-' upper(stimtype(D))];
    end
else
    IDstr = [icStr '-' iocStr '-' upper(stimtype(D))];
end

if isequal('full', flag),
    IDstr = [upper(name(D.ID.Experiment)) ', Rec ' idatStr '  <' IDstr '>'];
end

if isequal('status', flag),
    stat = status(D.ID.Experiment);
    icStr = num2str(stat.iCell);
    iocStr = num2str(stat.iRecOfCell+1);
    IDstr = [icStr '-' iocStr '-' upper(stimtype(D))];
end

