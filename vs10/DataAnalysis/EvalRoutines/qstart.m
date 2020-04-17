function qstart(varargin)
% personal MatLab startup file - by Abel
%
%%Pay attention to paths:
% datadir() should return the same as bdatadir() which is the path checked
% by databrowser. bdatadir is used in several dataset functions. Paths are
% saved in file: defaultDirs.SGSRsetup
%% -------------------- CHANGELOG --------------------
%VERSION: 0.01      Initial creation
%VERSION: 0.02      Adapted for C:\USR\UserName
%VERSION: 0.03      fix case sensitivity in mathlab 2014
%% -------------------- Configuration--------------------
% Users should only have to adapt the values in this "Configuration" field.

%My username:
% Previously, this global is set by <MatLabRoot>\local\qstart, i.e. the calling function.
% If empty take first argument as UserName.
global UserName				
if isempty(UserName)		
	UserName = varargin{1};
end
%Set local Disk
Disk = which('qstart');
%Early Configuration
EarlyPath = [Disk(1),':\EARLYrepo2015\'];
EarlyStim = [Disk(1),':\Early_StimDefLeuven2015\'];

%SGSR Configuration
SGSRPath = [Disk(1),':\SGSR\'];
StartUpDir = [SGSRPath, 'startupdir\'];
SetupInfo = [SGSRPath, 'SetupInfo\'];

WorkDir = [Disk(1),':\USR\', UserName];
mFilesDir = [WorkDir, '\mfiles\'];
localStartupDir = [WorkDir, '\MatlabStartup\'];
% DataDir = [WorkDir, '\rawData\'];
DataDir = [Disk(1),':\usr\adriaan\ExpData\'];
RapMcoDir = [Disk(1),':\Usr\', UserName, '\mco\'];
	
loadExperimentalTools = 0;
checkSGSRUpdates = 1;

%{
 Define incompatible paths
 - statslib: provides a "dataset" class which conflicts with the SGSR "dataset" class
 - Folders containing obsolete code && backupfiles
%}
Incompatible = { ...
    [matlabroot, '\toolbox\shared\statslib'],   ...
    [matlabroot, '\toolbox\dspblks'],...
    [SGSRPath, 'vs30\stimMENU\backup'],				...
	[SGSRPath, 'vs30\BramGeneral\_Obsolete'],		...
	[SGSRPath, 'vs30\BramGeneral\structtools\_obsolete'],	...
	[SGSRPath, 'vs30\Phoenix\Sys3\obsolete'],		...
    [SGSRPath, 'KQuest'], ...
    EarlyPath,...
    EarlyStim,...
    };

Paths = {       ...
    SGSRPath,   ...
    mFilesDir,  ...
    SGSRPath,	...
	RapMcoDir,	...
	localStartupDir,		...
    [SGSRPath, 'KQuest-Sandbox'] ...
    };
    
WarningsOff = {	...
	'MATLAB:dispatcher:InexactCaseMatch',	...
	'SGSR:Debug'
	};

	
	
%% -------------------- Main Function --------------------
% Users should not have to change this.

% supress warnings
warning('off', 'all')

% load the experimental tools?
file = [SetupInfo, 'localsysparam.SGSRsetup'];
if exist(file, 'file')
    conf = load(file, '-mat');
    conf.TDTpresent = loadExperimentalTools;
    save(file, '-struct', 'conf');
else
    TDTpresent = loadExperimentalTools;
    save(file, 'TDTpresent');
end

% Check for updates 
% if(checkSGSRUpdates)
%     matVersion = char(version);
%     cmd = ['!svn st -u --non-interactive --username updateCheck ',  ...
%             '--password updateCheck --no-auth-cache ',              ...
%             SGSRPath];
%     T=evalc(cmd);
%     if (matVersion(1) > 6)
%         k = regexp(T, '( \* )|(\n{0,1}! )','once');
%     else 
%         k = [strfind(T, ' * ') strfind(T, '!')];
%     end
%     if ~isempty(k)
%         disp('SGSR UPDATES AVAILABLE: Please update your SGSR version')
%     end
% end

% start SGSR
cd(StartUpDir);
save user.mat UserName; %Ugly hack to save username since javasetpath clears all vars 
startup;
load user.mat;			%Ugly hack to save username since javasetpath clears all vars 
delete user.mat;



% set current dir
if ~exist(mFilesDir,'dir')
    mkdir(mFilesDir);
end
cd(WorkDir);

%Add paths recursively
for n=1:length(Paths)
    addpath(genpath(Paths{n}));
end

% remove unwanted paths
for n=1:length(Incompatible)
    rmpath(genpath(Incompatible{n}));
end

% set the datadir 
dataDir(DataDir);

% express warnings
warning('on', 'all');

%remove unwanted warnings
for n=1:length(WarningsOff)
	warning('off', WarningsOff{n});
end



