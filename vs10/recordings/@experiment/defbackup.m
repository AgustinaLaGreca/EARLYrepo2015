function E=defbackup(E);
% experiment/defbackup - save a backup of experiment definition
%    defbackup(E) makes a backup of E by saving its definition to the
%    Foo.ExpDefBackup file, where Foo is name(E).


% load previous backups; append this one
FFN = [filename(E) '.ExpDefBackup']; 
if exist(FFN,'file'),
    BU = load(FFN,'-mat');
    BU = BU.BU;
else, % first backup
    BU = [];
end
BU = [BU E]; % append current one
save(FFN,'BU', '-mat');




