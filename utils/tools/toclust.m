function toclust(FN);
% toclust - copy EARLY mfile to CLUST

if ~isequal('KULAK', CompuName) && ~isequal('SIUT', CompuName) && ~isequal('Ee1285a', CompuName),
    error('ToClust only works at KULAK/SIUT/Ee1285a');
end
qq = which(FN,'-all');
if isempty(qq),
    error(['mfile ''' FN ''' not found.' ]);
end
if numel(qq)>1, error(['Multiple mfiles named ' FN '.']); end
FN = qq{1};
[DD FF] = fileparts(FN);
[DD1, DD2] = strtok(DD,'EARLY');
DestDir = fullfile('L:\', DD2);
if ~exist(DestDir,'dir'), 
    error(['Destination dir ''' DestDir ''' does not exist.']);
end
[SUCCESS,MESS,MESSAGEID] = copyfile(FN,DestDir);
error(MESS);


