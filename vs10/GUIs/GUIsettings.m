function GS = GUIsettings(FN);
% GUIsettings - settings for GUIs
%    GS = GUIsettings returns a struct whose fields contain settings for GUIs.
%    These settings are used by GUIpanel, ParamQuery, etc.
%  
%    GUIsettings('Foo') only returns the Foo category of the settings.

% only fill "by hand" when there is no persistent
persistent PGS; 
if isempty(PGS),
    PGS=local_spell_out; 
end

if nargin<1, FN=''; end

if isempty(FN),
    GS = PGS;
else,
    [FN, Mess] = keywordMatch(FN, fieldnames(PGS), 'GUIsetting category');
    error(Mess);
    GS = PGS.(FN);
end


%====================================
function GS = local_spell_out;
GS.general.Fontsize = 11; % fontsize for most stuff

% GUI figures
% figure properties
GS.GUIfigure.Color = [0.7 0.74 0.73]-0.05;
GS.GUIfigure.Units = 'normalized';
GS.GUIfigure.Position = [0.4 0.85 0.1 0.1];
GS.GUIfigure.KeyPressFcn = '1;'; % ignore keystrokes
GS.GUIfigure.IntegerHandle = 'on';
GS.GUIfigure.MenuBar = 'none';
GS.GUIfigure.Resize = 'off';
GS.GUIfigure.NumberTitle = 'off';
GS.GUIfigure.DockControls = 'off';
GS.GUIfigure.HandleVisibility = 'off';
% GUI panels
GS.GUIpanel.XYorigin = [5 23]; % position of Origin, serving as reference for first GUIpiece placed
GS.GUIpanel.XYlowrightMargin = [10 0]; % margins respected bij ADD method when updating GUIPanel extent
% uipanel properties
GS.GUIpanel.FontSize = GS.general.Fontsize;
GS.GUIpanel.FontWeight = 'bold';
GS.GUIpanel.BackgroundColor = GS.GUIfigure.Color;
GS.GUIpanel.ShadowColor = [0.45 0.55 0.45];

% ParamQuery objects
GS.ParamQuery.EditOkayColor = [1 1 1]; % background color of Edit uicontrols when okay
GS.ParamQuery.EditErrorColor = [1 0.7 0.7]; % background color of Edit uicontrols in case of error
GS.ParamQuery.ToggleErrorColor = [1 0.2 0.2]; % background color of pushbutton uicontrols in case of error
GS.ParamQuery.EditBlinkColor = 0.4*[0 1 1]; % idem blinking @ errors
GS.ParamQuery.EditBlinkDur = 0.05; % s blink duration
GS.ParamQuery.HeightFactor = 1.8; % blow-up factor relating heights of paramquery and its uicontrols
% ParamQuery uicontrol properties
GS.ParamQuery.Units = 'pixels';
GS.ParamQuery.FontSize = GS.general.Fontsize; 
GS.ParamQuery.FontName = get(0,'defaultuicontrolFontName');
GS.ParamQuery.FontUnits = 'points';
GS.ParamQuery.FontWeight = 'normal';

% ActionButton objects
GS.ActionButton.ExtentMatrix = [1 0 ; 0 1.5]; % matrix relating extent of string to extent of button
% ActionButton uicontrol properties
GS.ActionButton.Units = 'pixels';
GS.ActionButton.FontSize = GS.general.Fontsize; 
GS.ActionButton.FontName = get(0,'defaultuicontrolFontName');
GS.ActionButton.FontUnits = get(0,'defaultuicontrolFontUnits');
GS.ActionButton.BackGroundColor = 0.8*[1 1 0.95];

% Messenger objects
GS.Messenger.NeutralColor = [0 0 0]; % text color for neutral messages
GS.Messenger.WarnColor = 0.6*[1 0 1]; % idem warnings
GS.Messenger.ErrorColor = 0.6*[1 0 0]; % idem errors
GS.Messenger.BlinkColor = 0.4*[0 1 1]; % idem blinking @ errors
% Messenger uicontrol properties
GS.Messenger.FontSize = GS.general.Fontsize;
GS.Messenger.ForegroundColor = [0 0 0];

% AxesDisplay objects
% AxesDisplay axes properties
GS.AxesDisplay.Units = 'pixels';
GS.AxesDisplay.FontSize = GS.general.Fontsize;

% CycleList objects
% CycleList uimenu properties
GS.CycleList.ForegroundColor = [0 0 0];

% Pulldownmenu objects
% Pulldownmenu uimenu properties
GS.PulldownMenu.ForegroundColor = [0 0 0];





