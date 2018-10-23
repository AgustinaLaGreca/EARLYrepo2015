root = 'C:\EARLYrepo2015\GUIdefaults';
backupf = 'C:\EARLYrepo2015\GUIdefaults\0Examples\GUIdefaults';
cycleStr = 'CycleList*';
recentStr = 'most-recent*';

files = dir(root);
% which is a directory
dirIdx = [files.isdir];
subFolders = files(dirIdx);
subFolders=subFolders(~ismember({subFolders.name},{'.','..'}));

for f=1:length(subFolders)
folder = subFolders(f).name;
rootfolder = fullfile(backupf,folder);
subfiles = dir(rootfolder);
subfiles = {subfiles.name};
% no all folders have both files
try
cycle = subfiles{~cellfun(@(x) isempty(regexp(x,cycleStr,'once')),subfiles)};
recent = subfiles{~cellfun(@(x) isempty(regexp(x,recentStr,'once')),subfiles)};
try
    origin = {fullfile(rootfolder,cycle),fullfile(rootfolder,recent)};
    dest = {fullfile(root,folder,cycle),fullfile(root,folder,recent)};
    cellfun(@copyfile,origin,dest);
    disp(['copied ' folder]);
catch
    disp([folder ' error copying in this folder'])
end
catch
    disp(['Empty folder ' folder])
end
end

