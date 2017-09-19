function versiondir = EarlyUpgrade(versiondir);
% EarlyUpgrade - upgrade EARLY files from zip file in startup dir
%
%    EarlyUpgrade, which is called at startup time via Startup, checks
%    whether there is a zipped version of Early in the startup directory.
%    If so, the present Eraly installation is backed up, and the Winzip
%    file is unpacked to produce the current Early installation.
%
%    See also startup, earlyver.

Mess = ''; % default: no problems
EarlyRoot = fileparts(versiondir); % parent dir

% -----find out if an update is due
zipfilelist = dir('*.zip');
if length(zipfilelist)>1,
    warning('Multiple ZIP file in startup dir. Clean up to restore proper processing of upgrades.');
elseif isempty(zipfilelist), % no upgrade
    return;
elseif isempty(strfind(zipfilelist(1).name,'EARLY')),
    warning(['Zipfile ''' zipfilelist(1).name ''' in startup directory ignored - apparently not an EARLY upgrade.' ]);
    return;
end

% -----backup current versiondir by renaming it
[PP NN] = fileparts(versiondir);
vsBackup = fullfile(PP, ['BACKUP_of_VSdir']);
if exist(vsBackup, 'dir'), % remove any previous backups of versiondir
    [Success,Mess] = rmdir(vsBackup, 's');
    if ~Success, error(Mess); end
end
cmd = ['!move "' versiondir '" "' vsBackup '"'];
eval(cmd);
if exist(versiondir,'dir') || ~exist(vsBackup,'dir'),
    error(['Renaming versiondir ''' versiondir ''' failed.']);
end
% store name of orginal vsdir (reason for generic naming of backupdir is to allow zipwork to avoid including it in backups)
[fid, Mess] = fopen(fullfile(vsBackup,'README.TXT'), 'at');
error(Mess);
fprintf(fid, '%s BACKUP of %s', datestr(now), versiondir);
fclose(fid);

% ----unzip the zip file
ZipfileName = fullfile(pwd, zipfilelist.name);
exe = '"C:\Program Files\WinZip\wzunzip" -don '; % d=use dirs; o=overwrite mode; n=only newer files (to prevent overwriting setup files!)
cmd = [exe ZipfileName ' ' EarlyRoot];
disp('        .....updating EARLY ..... please wait .... ');
[status, result] = dos(cmd);
if status, % winzip returns ~0 if an error occurred. Try to recover
    try,
        disp('Attempting to recover from failed update ..')
        [Success,Mess] = rmdir(versiondir, 's');
        cmd = ['!move "' vsBackup '" "' versiondir '"'];
        eval(cmd);
        h = errordlg({'Failed attempt to updata EARLY toolbox' 'Check Zipfile in EARLY\startupdir!'}, 'No update for you', 'modal');
        waitfor(h);
    catch,
        h = errordlg({'Failed recovery from failed attempt to updata EARLY toolbox' 'EARLY initialization will be interrupted.'}, 'Now it''s a real mess!', 'modal');
        waitfor(h);
        clc; disp(blanks(4)');
        error('Unable to recover from failed EARLY upgrade. Restore backup of versiondir, remove zipfile from startupdir, and restart EARLY.');
    end
    return;
end
% ----versiondir may have changed due to upgrade; check it again
versiondir = findversiondir(EarlyRoot); 
%clc;
disp(' ')
disp('---EARLY updated---');
disp(' ')

% Remove the zipfile, after logging its unzipping
UpgradeLog = fullfile(prefdir(1), 'EARLY_update.log');
[fid, Mess] = fopen(UpgradeLog, 'at');
error(Mess);
fprintf(fid, '%s Unzipped %s\n', datestr(now), zipfilelist.name);
fclose(fid);
delete(ZipfileName);






