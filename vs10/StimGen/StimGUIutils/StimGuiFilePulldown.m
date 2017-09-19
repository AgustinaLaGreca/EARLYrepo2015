function Pfile=StimGuiFilePulldown(src, dum, ID, varargin);
% StimGuiFilePulldown - creator & callback File Pulldown menu of StimGUI
%    P=StimGuiFilePulldown returns the PulldownMenu object to be
%    included in a StimGUI. It contains the standard repertory of menu
%    items: 
%      Open file and read parameters 
%      Save parameters in file
%      Empty cycle List
%      CycleList of Recently read parameters files
%
%    StimGuiFilePulldown also serves as the callback function of the very
%    menu items created by it.
%
%    See also StimGUI, StimGuiEditPulldown, StimGuiViewPulldown.

if nargin<1, % create call
    callme = fhandle(mfilename); % function handle to this function 
    Pfile=pulldownmenu('File','&File');
    Pfile=additem(Pfile,'&Save parameters', {callme 'save'}, 'accelerator', 'S');
    Pfile=additem(Pfile,'&Open parameter file', {callme 'read'}, 'accelerator', 'O');
    Pfile=additem(Pfile,'&Empty File List', {callme 'emptylist'});
    CL=cycleList('StimparamDefaults', 10, callme); % clicking an item calls this fcn with 4 input args; see CycleList
    Pfile=additem(Pfile,CL);
else, % callback
    figh = parentfigh(src); % GUI figure handle
    if isequal('save', ID), 
        % ==========prompt user for param file & save=========
        GUIgrab(figh,'?',1); % last arg: 1 = do display GUImessage
    elseif isequal('read', ID), 
        % =========prompt user for param file & read=========
        [Mess, Ignored, FileName]=GUIfill(figh,'?',1); % last arg: 1 = do display GUImessage
        if ~isempty(FileName), % add filename to cycleList
            [PP FileName ee] = fileparts(FileName); % strip off dir & extension
            % first retrieve cycle list C, then add 
            hP = get(src, 'parent'); % parent of item is uimenu
            P = get(hP,'userdata'); % Pulldown object
            C = getCycleList(P,'StimparamDefaults');
            C = additem(C,FileName);
        end
    elseif isequal('emptylist', ID), 
        % ============clear the cycle list============================
        % first retrieve cycle list C
        hP = get(src, 'parent'); % parent of item is uimenu 
        P = get(hP,'userdata'); % Pulldown object
        iC = find(P.CycleItems, 'StimparamDefaults');
        C = P.CycleItems(iC(1));
        % clear its contents
        C = rmitem(C);
    elseif isCycleList(ID), %4th arg passed to callback is item (see cycleList/select)
        %=============read params from file of cycleList item =====================
        C = ID; % the up-to-date cycleList
        item = varargin{1}; FileName = item.Label;
        [Mess, Ignored, FileName]=GUIfill(figh,FileName,1);
        if isempty(FileName), % remove item from list as it appears obsolete
            C=rmitem(C,item.Label);
        end
    else, 
        error('Invalid callback.');
    end
end
