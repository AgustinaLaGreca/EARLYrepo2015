function [EarlyRoot, versiondir] = startup(varargin)
% startup - initialization for EARLY toolbox
%   startup is called at MatLab initialization time by starting MatLab in
%   the directory that contains the startup.m file.
%   Alternatively, you may call startup from your own initialization script
%   or function.
%
%   Startup finds out what the current directory is and sets the EARLY path
%   accordingly. It also sets some other default Matlab settings.
%
%   [EarlyRoot, versiondir] = startup returns the root dir and version dir
%   of the EARLY toolbox.
%
%   Note: startup first restores Matlab's default path and then adds the
%   EARLY path dirs. Any user-defined paths should thus be added *after*
%   calling startup.
%
%   See also versiondir.


warning('off','MATLAB:dispatcher:InexactMatch');
warning('off','MATLAB:dispatcher:InexactCaseMatch');
set(0,'FormatSpacing','compact');
set(0,'DefaultFigureCreateFcn',@CreateFigCallback); % enable line pasting

% restore Matlab path at installation and add current dir
restoredefaultpath;
rmpath(fullfile(matlabroot, 'toolbox', 'shared', 'statslib')); % removes dataset.m constructor in shared/statslib
rmpath(fullfile(matlabroot, 'toolbox', 'nnet', 'nnet','nndatafun')); % removes minmax.m in nnet

% find out EARLY's root dir; by convention, startup.m (this file) resides
% in EarlyRoot\startup
startupdir = fileparts(mfilename('fullpath'));
addpath(startupdir); % so that startup.m is in the path
EarlyRoot = fileparts(startupdir); % parent dir

% find cur version dir vsxx: % select dirs named vsxx and find highest xx
versiondir = findversiondir(EarlyRoot); % findversiondir must reide in startupdir; we have no EARLY path yet

% update to new version if present.
versiondir = EarlyUpgrade(versiondir);

% not all dirs under versiondir need to be in path; also some path settings
% may be conditional on local setup features, users, etc. This is handled
% by EarlyPath in the init/ dir.
cd(fullfile(versiondir, 'init'));
EarlyPath(versiondir);

% kick random generators
SetRandState;

% set cur dir to harmless dir outside path
cd(fullfile(EarlyRoot));

if nargin <1
    username(100); % Ensure a username is known: Username will prompt the user if not.
end

disp('----EARLY toolbox initialized----'); disp(' ');

%by abel print message of the day
motdTxt = fullfile(EarlyRoot, 'doc', 'motd.m');
if exist(motdTxt, 'file')
    disp('<---------------------- MOTD ---------------------->'); 
    type(motdTxt)
    disp('<---------------------- %%%% ---------------------->'); 
end

if nargout<1, clear EarlyRoot; end % turn off echoing @ startup

