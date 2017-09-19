function S=GUIgrab(figh,FN,doDisp);
% GUIgrab - get strings from GUI edits and toggles without interpreting them
%    S=GUIgrab(h) returns all strings of edit controls and toggle buttons
%    in struct S. Unlike GUIval, GUIgrab returns the string values as they
%    are currently displayed on the GUI, without converting them to
%    numbers, checking syntax, or validity of the values. Empty edit 
%    controls are left out of S, as are disabled toggle buttons (toggles 
%    are enabled/disabled by right-clicking). GUIgrab is used to store a
%    set of GUI values for future use (Default settings). The values
%    collected in S can be re-displayed in a GUI by GUIfill.
%
%    S=GUIgrab(h, 'Foo') saves the grabbed param info to GUIdefaults file
%    Foo. For info on directories, see GUIfill.
%
%    S=GUIgrab(h, '?') prompts for a filename and saves the grabbed param
%    to the specified GUIdefaults file.
%
%    S=GUIgrab(h, '>') appends the the grabbed param info to history
%    file of this GUI type. See GUIfill for details.
%
%    See also GUIval, GUIfill, paramquery/grab.

Nhistory = 20; % max number of stored histories
if nargin<2, FN=[]; end
if nargin<3, doDisp=0; end % default: don't display message

figure(figh);
S.GUIname = getGUIdata(figh,'GUIname');
S.EnableStates_ = [];
%GUImessage(figh,'','neutral'); % clear any previous messages
Q = getGUIdata(figh,'Query', []); % array of queries.
% read queries one by one. 
for ii=1:numel(Q),
    q = Q(ii);
    [EditStr, UnitStr, enab]=grab(q);
    if ~isempty(EditStr), 
        S.([q.Name 'Edit']) = EditStr;
        S.EnableStates_.(q.Name) = enab;
    end
    if ~isempty(UnitStr), 
        S.([q.Name 'Unit']) = UnitStr;
    end
end

if ~isempty(FN) && ~ischar(FN),
    error('FN input arg must be filename or ''>''.')
elseif ~isempty(FN) && ~isequal('>', FN), % ----------save to guidef file
    DD = GUIdefaultsDir(S.GUIname, 'create'); % GUIdef dir of this GUI type; make sure it exists
    if isequal('?', FN), % prompt for file
        fn = uiputfile(fullfile(DD,'*.GUIdef'), ['Saving to file:  specify filename for ' S.GUIname ' parameters']);
        if isequal(0,fn), return; end % user canceled
        [dum, FN, dum2] = fileparts(fn); % strip off dir & ext
    end
    FN = FullFilename(FN, DD, '.GUIdef'); % full file name, including extension
    grabStruct = S;
    save(FN, 'grabStruct', '-mat');
    [DD, fn, EE] = fileparts(FN);
    if doDisp, GUImessage(figh, ['Parameters saved to file ''' fn '''']); end
elseif isequal('>', FN), % %-------------save to history file
    DD = GUIdefaultsDir(S.GUIname, 'create'); % GUIdef dir of this GUI type; make sure it exists
    FN = FullFilename('most-recent', DD, '.GUIhistory');
    if exist(FN,'file'), % append
        S_old = load(FN,'-mat'); S_old = S_old.grabStruct;
        % only store if different from last one
        if ~isequal(S_old(end), S),
            grabStruct = [S_old S];
        end
        if length(grabStruct)>Nhistory, % limit the length of the stored history
            grabStruct = grabStruct(end-Nhistory+1:end);
        end
    else, % create
        grabStruct = {S};
    end
    save(FN, 'grabStruct', '-mat');
end





