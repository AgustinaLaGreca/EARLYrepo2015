function gomaster()
%% GOMASTER - Checkout EARLYrepositories to branch master for experiments
% It should be executed before starting a set of experiments in order to
% ensure that EARLY is on the most recent stable version, which is the last
% commit of branch master
% It commits the current state of the develop branch if neccesary and then
% checkout
% Located in startupdir in order to be added to the path at starting up

% These var contain the bat files that perform the desired actions on git
% They have absolute paths to repositories
commit = 'tempCommit.bat';
checkout = 'checkoutmaster.bat';
pull = 'pull.bat';

startupdir = fileparts(mfilename('fullpath'));
EarlyRoot = fileparts(startupdir); % parent dir
stimRoot = fileparts(EarlyRoot);
stimRoot = fullfile(stimRoot,'EARLY_StimDefLeuven2015');

batPath = fullfile(EarlyRoot,'gitUtils');
commFilename = fullfile(batPath,commit);
checkFilename = fullfile(batPath,checkout);
pullFilename = fullfile(batPath,pull);

mbranch = 'master';

% Check in which branch EARLY is
cmd1 = ['cd ' EarlyRoot];
cmd1b = ['cd ' stimRoot];
cmd2 = 'git rev-parse --abbrev-ref HEAD';
cmd = [cmd1 ' & ' cmd2];
cmdb = [cmd1b ' & ' cmd2];
[~, inbranch1] = system(cmd);
[~, inbranch2] = system(cmdb);

if ~strcmp(inbranch1(1:end-1),mbranch) || ~strcmp(inbranch2(1:end-1),mbranch)
    [~, ~] = system(commFilename);
    [~, ~] = system(checkFilename);
    [~, pullmsg] = system(pullFilename);
    [~, inbranch1] = system(cmd);
    [~, inbranch2] = system(cmdb);
    bnchmsg = ['You are now in ' inbranch1(1:end-1) ' and ' inbranch2(1:end-1) ' branch.'];
    msg = {bnchmsg,pullmsg};
    fprintf('%s\n',msg{:});
else
    disp('Already in master branch.');
end
