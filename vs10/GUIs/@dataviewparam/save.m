function [okay, Filename] = save(P, Filename, DefName);
% dataviewparam/save - save dataviewparam object to disk
%   save(P, 'Foo') saves P to file Foo.dvparam in the 'Dataviewparam/XXX'
%   subdirectory of the GUIdefaultsDir() directory, where XXX is the name
%   of P's dataviewer. When trying to overwrite the file named def.dvparam,
%   the user is asked to confirm his intention to change the default
%   dataviewparam object for XXX.
%
%   [Saved, FN] = save(P, '?') prompts the user for a filename and 
%   then saves, unless the user cancels. Saved will be True if the P was 
%   saved, False,  otherwise. The full filename is returned in FN (empty
%   when not saved).
%   
%   save(P, '?', 'Foo') prompts the user for a filename while using a
%   default name Foo, then saves unless the user cancels.
%   
%   See also GUIdefaultsDir, dataviewparam/savedir, dataviewparam/load.

if nargin<3, DefName=''; end
okay = 0;
if isvoid(P),
    error('Cannot save a void dataviewparam object.');
end

doSetDefault = 0;
sdir = savedir(P);
while isequal('?', Filename), %prompt
    Filefilter = fullfile(sdir, '*.dvparam');
    if ~isempty(DefName), DefName = fullfile(sdir, DefName); end
    Filename = uiputfile(Filefilter, 'Save to where?', DefName);
    if isequal(0, Filename), Filename=''; 
        return; % user cancelled -> quit
    end
    Filename = strtok(Filename, '.'); % remove extension
    % if the user selected 'def', ask for confirmation
    if isequal('def', lower(Filename)),
        DV = char(dataviewer(P));
        QQQ= questdlg({['Your choice will change the default settings for ' DV ' parameters'], ...
            'Are you sure you want to change this default?'}, ...
            ['Change ' DV ' defaults??'], ...
            'Yes', 'No', 'No');
        if isequal('No', QQQ), Filename = '?'; % this keeps the prompting going
        else, doSetDefault = 1;
        end
    end
end

XFN = [Filename, '.dvparam']; % extended filename; note that extensions includes viewer
CFN = fullfile(sdir, XFN);
%CFN
CacheParam = 1; % not a true cache, just save once
putcache(CFN, 1, CacheParam, P);
if doSetDefault, setdefault(P); end
okay = 1;
Filename = CFN; % return full file name as described in help text.



