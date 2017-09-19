% quick copyscript
DestDir = 'G:\Thomas\data\sortTKdata\tempABFstore';
qq = load(fullfile(DestDir, 'tempABFdirList.mat'), '-mat');
Copied = qq.tempABFdir(3:end); % remove .\ and ..\
CopiedNames = {Copied.name};
Pool = dir('D:\USR\Marcel\sortTKdata\tempABFstore\*.abf');
ineeded = [];
for ii=1:numel(Pool),
    p = Pool(ii);
    if isempty(strmatch(p.name, CopiedNames)),
        ineeded = [ineeded, ii];
    end
end
%ineeded
ToBeCopied = Pool(ineeded);
for ii=1:numel(ToBeCopied),
    tb = ToBeCopied(ii);
    display(tb.name);
    copyfile(fullfile('D:\USR\Marcel\sortTKdata\tempABFstore',tb.name), DestDir);
end



