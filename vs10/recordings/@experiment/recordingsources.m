function [Src, Nm] = recordingsources(E, OnlyActiveOnes, dashboardGUIval);
% experiment/recordingsources - recording sources defined for experiment
%   [Src, Nm] = recordingsources(E) returns struct Src containing the
%   info on the recording source defined for experiment E. Nm is a cell
%   array of their names, i.e., Nm = fieldnames(Src);
%   If E is void, Src is a void struct and Nm equals {}.
%
%   Recordingsources(E, 1) returns only the "active" sources, i.e.,
%   those selected in the edit GUI of the experiment, where inactive 
%   channels are specfied by buttons showing '-'.
%   The default is to show all sources, regardless of their being active.
%
%   Recordingsources(E, 1, RS) further filters the sources using the struct
%   RS resulting from GUIval-ling the dashboard. In the dashboard one can
%   switch on and off the recording of each channel. Disabling a channel is
%   done by tggling to a '-' value.
%
%   See Experiment/edit, Experiment/isvoid.

[OnlyActiveOnes, dashboardGUIval] = arginDefaults('OnlyActiveOnes/dashboardGUIval', 0,[]); % def: all sources; no dashboard selection

if isvoid(E),
    Src = struct();
else,
    Src = E.Recording.Source;
    if isempty(Src),
        Src = struct();
    end
end
Nm = fieldnames(Src); Nm = Nm(:).';


if OnlyActiveOnes, % filter while using the DataType property
    for ii=1:numel(Nm),
        fn = Nm{ii};
        if isequal('-', Src.(fn).DataType), % eliminate it
            Src = rmfield(Src, fn);
        end
    end
    Nm = fieldnames(Src); Nm = Nm(:).';
end

if ~isempty(dashboardGUIval),
    for ii=1:numel(Nm),
        fn = Nm{ii};
        ToggleValue = trimspace(dashboardGUIval.(fn));
        if isequal('-', ToggleValue), % eliminate it
            Src = rmfield(Src, fn);
        end
    end
    Nm = fieldnames(Src); Nm = Nm(:).';
end
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    




