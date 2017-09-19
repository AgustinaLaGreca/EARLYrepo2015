function StimGuimode(figh,Mode,gmess, gmessmode);
% StimGuimode - set enable/disable mode of action buttons of stimulus GUI
%    StimGuimode(figh,Mode) sets the button mode to Mode, which is one of
%    'Ready' 'Busy' 'Check' 'Play' 'PlayRecord' 'Stop'. 
%    (Abbreviations & case mismatches are allowed, see keywordMatch.)
%    This affects the enable mode of the buttons, their hightlighting, and 
%    the ability to close the GUI figure.
%
%    StimGuimode(figh,Mode, Mess, MessMode) also displays text message Mess
%    using GUImessage. The message mode MessMode (see GUImessage) defaults
%    to 'neutral'.
% 
%    See also StimGUIaction.
if nargin<3, gmess = nan; end % indicates absence of message ('' would be bad choice ;)
if nargin<4, gmessmode = 'neutral'; end
[Mode, Mess] = keywordMatch(Mode,{'Check' 'Busy' 'Ready' 'Play' 'PlayRecord' 'Stop'}, 'Mode argument');
error(Mess);

A = getGUIdata(figh, 'ActionButton');
highlight(A,'default'); % set buttons to their default colors
switch Mode,
    case {'Check', 'Busy', 'Stop'}, % disable all buttons & prevent closing the figure
        enable(A,0);
        GUIclosable(figh,0); % not okay to close GUI
        % color Check or Stop buttons
        if isequal('Check', Mode),
            highlight(A('Check'),[0 0 0.7]);
        elseif isequal('Stop', Mode),
            highlight(A('Stop'),[0.5 0.15 0]);
        end
    case 'Ready', % enable all buttons except Stop; 
        enable(A,1);             
        enable(A('Stop'),0);
        % disable PlayRec if no experiment is ongoing
        enable(A('PlayRec'), ~isvoid(experiment('current')));
        % only enable Play if D/A is possible
        enable(A('Play', 'PlayRec'), CanPlayStim);
        GUIclosable(figh,1); % okay to close GUI
    case {'Play' 'PlayRecord'}, % disable all buttons except Stop
        enable(A,0);
        enable(A('Stop'),1);
        GUIclosable(figh,0); % not okay to close GUI
        % color Play or PlayRec buttons
        if isequal('Play', Mode),
            highlight(A('Play'),[0 0.7 0]);
        elseif isequal('PlayRecord', Mode),
            highlight(A('PlayRec'),[0.85 0 0]);
        end
end
% display GUI message, if any.
if ~isnan(gmess),
    GUImessage(figh,gmess,gmessmode);
end

