function GUIreaper(figh,flag);
% GUIreaper - access or create GUI reaper
%    GUIreaper(figh, 'create') creates an invisible UIcontrol element that
%    serves to signal functions to close the figure. There can only exist
%    one reaper per figure. 
%
%    GUIreaper(figh, 'wait') calls waitfor on the figure's reaper, i.e., it
%    blocks execution until a GUIreaper(figh,'reap') has been called. If no
%    reaper has been defined yet, GUIreaper(figh, 'create') iscalled
%    implicitly.
%    
%    GUIreaper(figh, 'reap') triggers blocked execution by deleting the
%    hidden reaper uicontrol.
%
%    See also GUIclose, gcg.

if ~istypedhandle(figh,'figure'), error('Handle does not belong to a figure.'); end

[flag, Mess] = keywordMatch(lower(flag), {'create' 'wait' 'reap'},'flag');
error(Mess);

switch flag
    case 'create',
        hr = getGUIdata(figh, 'reaper',[]);
        if ~isempty(hr),
            error('Cannot create multiple reapers for a single figure.');
        end
        hr = uicontrol(figh, 'tag', 'reaper', 'visible','off');
        setGUIdata(figh, 'reaper', hr);
    case 'wait',
        hr = getGUIdata(figh, 'reaper',[]);
        if isempty(hr), % no reaper exists; create one
            GUIreaper(figh,'create');
            hr = getGUIdata(figh, 'reaper',[]);
        end
        waitfor(hr);
    case 'reap',
        hr = getGUIdata(figh, 'reaper',[]);
        if isempty(hr),
            error('No reaper has been created for this figure; cannot reap.');
        end
        delete(hr);
end




