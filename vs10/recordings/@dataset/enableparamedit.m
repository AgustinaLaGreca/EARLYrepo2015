function figh = enableparamedit(DS, P, figh);
% Dataset/enableparamedit - enable editing of dataviewer parameters
%    figh = enableparamedit(DS, P, figh) provides uimenu item and a toolbar 
%    pushbutton to the figure with handle figh. The menu item and the 
%    toolbar pushbutton both launch the editing of the
%    parameters of the dataviewer in dataviewparam object P. After 
%    confirming the edit, the figure is the analysis is redone using the 
%    new parameter values, and the figure is refreshed.
%
%    See dataviewparam, dataset/listdataviewer, dataset/isdataviewer, 
%    dataset/dotraster.

if ~isa(P, 'dataviewparam'),
    error('Second input argument must be a Dataviewparam object.')
end

if ~isSingleHandle(figh, 'figure'),
    error('Figh argument must be handle to existing figure.');
end
% don't enable twice
ie = getGUIdata(figh, 'paramEditEnabled', false); 
if ie, return; end % already enabled - done

% store DS and P in GUIdata
setGUIdata(figh, 'Dataset', DS);
setGUIdata(figh, 'DataviewParam', P);
setGUIdata(figh, 'paramEditEnabled', true); % notify that param editing is enabled
% add uimenu & item
hParam = uimenu(figh, 'label', '&Params', 'Tag', 'ParameterMenu', 'HandleVisibility', 'off');
hEditParam = uimenu(hParam, 'label', 'Edit', 'Accelerator', 'Q', ...
    'callback', @(Src, Ev)local_edit(figh), 'Tag', 'EditParamsMenuItem', ...
    'HandleVisibility', 'off');

% =========================================
function local_edit(figh);
% edit the parameters and refresh view
P = getGUIdata(figh, 'DataviewParam');
Pe = edit(P);
if isvoid(Pe), % user cancelled
    return;
end
DS = getGUIdata(figh, 'Dataset');
dataview(DS, Pe, figh);
setGUIdata(figh, 'DataviewParam', Pe);






