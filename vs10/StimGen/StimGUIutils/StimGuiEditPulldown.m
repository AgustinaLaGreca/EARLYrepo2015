function P=StimGuiEditPulldown(src, dum, ID, varargin);
% StimGuiEditPulldown - creator & callback Edit Pulldown menu of StimGUI
%    P=StimGuiEditPulldown returns the PulldownMenu object to be
%    included in a StimGUI. It contains the standard repertory of menu
%    items: 
%      Empty all edit fields
%      Undo changes since last Check & Update
%      Redo undone change
%
%    StimGuiEditPulldown also serves as the callback function of the very
%    menu items created by it.
%
%    See also StimGUI, StimGuiFilePulldown, StimGuiViewPulldown.

if nargin<1, % create call
    callme = fhandle(mfilename); % function handle to this function 
    P=pulldownmenu('Edit','&Edit');
    P=additem(P,'&Empty all fields', {callme 'empty'}, 'accelerator', 'E');
    P=additem(P,'&Undo changes since last Check', {callme 'undo'}, 'accelerator', 'Z', 'userdata',0, 'tag','Undoer');
    P=additem(P,'&Redo undone changes', {callme 'redo'}, 'accelerator', 'Y');
    
else, % callback
    figh = parentfigh(src); % GUI figure handle
    if isequal('empty', ID), 
        % ==========empty all edit fields=========
        Q=getGUIdata(figh, 'Query');
        clearedit(Q);
    elseif isequal('undo', ID), 
        % =========undo changes since last Check=========
        HC= get(src,'userdata');
        if isequal(0,HC), % initialize
            HistCount=0; 
            HistDate=localHistoryDate(figh);
            HC = CollectInStruct(HistCount,HistDate);
            set(src,'userdata',HC);
        end
        % has history been updated since last undo call? If so, reset counter
        HistDate=localHistoryDate(figh);
        if ~isequal(HC.HistDate, HistDate), % reset
            HC.HistCount=1;
            HC.HistDate=HistDate;
        end
        HC.newestState = GUIgrab(figh);  % store current GUI state
        GUIfill(figh,HC.HistCount); 
        HC.HistCount=HC.HistCount+1;
        set(src,'userdata', HC);
    elseif isequal('redo', ID), 
        hu = findobj(figh, 'tag','Undoer');
        HC= get(hu,'userdata');
        if isequal(0,HC), return; end % no recorded history
        HC.HistCount = HC.HistCount-1;
        if HC.HistCount<0,HC.HistCount=0; end
        if HC.HistCount==0, GUIfill(figh,HC.newestState);
        else, GUIfill(figh,HC.HistCount); 
        end
        set(hu,'userdata', HC);
    else, 
        error('Invalid callback.');
    end
end
%========================
function hd=localHistoryDate(figh);
% date number of hostory file for this GUIname
GUIname = getGUIdata(figh,'GUIname');
DD = GUIdefaultsDir(GUIname, 'create'); % GUIdef dir of this GUI type; make sure it exists
FN = FullFilename('most-recent', DD, '.GUIhistory');
if exist(FN,'file'),
    qq=dir(FN);
    hd = qq(1).datenum;
else,hd=0;
end


