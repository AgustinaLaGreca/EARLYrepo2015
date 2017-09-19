function F = folder(E);
% experiment/folder - data folder (directory) of experiment on this computer
%    folder(E) returns the location of experimental files of experiment E
%    on the current computer. It is located by recursively searching for a
%    directory named XXX, where XXX=name(E), in the subfolders of
%    parentfolder(E). The first match is returned.
%     
%    See also experiment/parentfolder, experiment/filename, experiment/save.

persistent FFFF

if isvoid(E),
    error('Void experiments do not have a folder.');
end

% local cache to avoid having to search subfolders all the time
if ~isempty(FFFF),
    F = getFieldOrDefault(FFFF, name(E), pi); % impossible default value
    if ~isequal(pi, F), return; end
end
% real cache
[F, CFN, CP] = getcache(mfilename, {parentfolder(E) name(E)});
if isempty(F),
    % experiment directory is of experiment FOO is sub-(sub)-dir of
    % parentfolder that contains FOO.ExpDef file.
    %F = searchInPath(name(E), genpath(parentfolder(E)), 'dir');
    F = local_search(name(E), parentfolder(E));
    % store in cache if not empty; don't store empty - new experiments may be
    % defined in the mean time.
end

if ~isempty(F), % caches 
    putcache(CFN, 1e3, CP, F);
    FFFF.(name(E)) = F;
end

%======================================================
function F = local_search(Name, Pdir);
F = [];
Tg = fullfile(Pdir, Name);
if exist(Tg, 'dir'),
    F = Tg;
    return;
end
Edir = dir(Pdir);
Edir = Edir(3:end); % exclude .\ and ..\
Edir = Edir([Edir.isdir]);
for ii=1:numel(Edir),
    Tg = fullfile(Pdir, Edir(ii).name, Name);
    if exist(Tg, 'dir'),
        F = Tg;
        return;
    end
end










