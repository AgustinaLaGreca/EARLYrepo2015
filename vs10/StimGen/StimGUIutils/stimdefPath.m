function P = stimdefPath(kw, varargin);
% stimdefPath - list of folders containing stimulus definitions
%     stimdefPath returns a semicolon-separated path for locating stimulus
%     definitions. By convention, the stimulus definition for 
%     stimulus Foo is in a file names stimdefFoo in this path. 
%     The first folder of stimdefPath is always Early's native stimdef 
%     folder, i.e., stimdefDir. The factory setting (see below) has 
%     stimdefPath equal to stimdefDir, without any additional folders.
%    
%     stimdefPath('append', 'Dir1', 'Dir2', ...) appends additional
%     directories to the factory settings. These persist across 
%     Matlab sessions. Full directory names must be used; specify all of 
%     them in a single call, because each subsequent call disables any
%     previously added dirs. Note that these directories will be added to 
%     the Matlab path at a position afther the Eary-specific path and
%     before the standerd Matlab path.
%
%     stimdefPath('factory') resets the stimdef path to stimdefDir only and
%     removes any previously specified stimdef directories from the Matlab
%     path.
%
%     See also EarlyPath, StimCollection, stimdefDir, stimdefFS.

%<------------------ CHANGE LOG ------------------>
% by: Abel on: di apr  5 13:57:45 CEST 2016
% Replaced hard coded ';' by filesep

persistent PP

kw = arginDefaults('kw', 'get');
SFN = 'LocalStimuli'; % setup filename

switch lower(kw),
    case 'get',
        if isempty(PP),
            PP = FromSetupFile(SFN, 'StimdefPathExtension', '-default', '');
        end
    case 'append',
        % undo previous dirs
        stimdefPath factory;
        % insert dirs into path (will fail if any dir is nonexistent)
        EarlyPath('insert', varargin{:}); % pre-matlab, post-early
        StimdefPathExtension = cell2words(varargin, pathsep);
        ToSetupFile(SFN, StimdefPathExtension);
        PP = FromSetupFile(SFN, 'StimdefPathExtension', '-default', '');
    case 'factory',
        % remove current list from Matlab path
        qq = FromSetupFile(SFN, 'StimdefPathExtension', '-default', '');
        qq = Words2cell(qq,pathsep);
        for ii=1:numel(qq), rmpath(qq{ii}); end
        % replace stored path extension by ''
        StimdefPathExtension = '';
        ToSetupFile(SFN, StimdefPathExtension);
        PP = '';
    otherwise,
        error('Invalid keyword.');
end
%set local disk
Disk = which('startup');
if ~strcmp(Disk(1),'C'),PP=strrep(PP,'C:\',Disk(1:3));end
P = [stimDefDir pathsep PP];



