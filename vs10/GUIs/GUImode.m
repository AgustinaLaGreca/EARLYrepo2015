function GUImode(figh, Mode,gmess,gmessmode)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
if nargin<3, gmess = inf; end % indicates absence of message - '' would be bad choice.
if nargin<4, gmessmode = 'neutral'; end
[Mode, Mess] = keywordMatch(Mode,{'Busy' 'Ready' 'Play' 'PlayRecord' 'Replay' 'Stop'}, 'Mode argument');
error(Mess);
A = getGUIdata(figh, 'ActionButton');
A_ExpStat = A('NewUnit', 'NewElectrode', 'Note', 'FinishExp', 'EditExp'); % buttons changing  the Experiment status
Q = getGUIdata(figh, 'Query', ParamQuery()); if numel(Q) == 1 && isvoid(Q), Q = Q([]); end
Exp = current(experiment);
[dum RecSrc] = recordingsources(Exp);
Qmeasure = Q(RecSrc{:}); % queries having to do w recordings
hasExp = ~isvoid(Exp); % true if experiment is going on
if isvoid(Exp), ExpStr = ' (no experiment)'; else, ExpStr = [' -- Experiment ' name(Exp) ' by ' experimenter(Exp)]; end
set(figh, 'name', ['Dashboard' ExpStr]);

hasStim = isSingleHandle(getGUIdata(figh, 'StimGUIhandle'));
switch Mode,
    case {'Busy', 'Stop'}, % disable all buttons & prevent closing the figure
        enable(A,0); enable(Q,0);
        GUIclosable(figh,0); % not okay to close GUI
        % color Check or Stop buttons
        if isequal('Stop', Mode),
            highlight(A('Stop'),[0.5 0.15 0]);
        end
    case 'Ready', % enable all buttons except Stop; okay to close GUI; recording queries depend on experiment status
        enable(A,1);  enable(Q,1); enable(Qmeasure, hasExp);
        %enable(A('StimSpec'), hasExp);
        enable(A('Stop'),0);
        % only enable Play if D/A is possible; only enable PlayRec when an experiment is ongoing
        enable(A('Play'), CanPlayStim && hasStim);
        enable(A('PlayRec'), CanPlayStim && hasStim && canrecord(Exp));
        enable(A_ExpStat, hasExp); % Exp status may only be changed when an Exp has been defined
        enable(A('NewExp', 'ResumeExp'), ~hasExp);
        highlight(A,'default');
        GUIclosable(figh,1); % okay to close GUI
    case {'Play' 'PlayRecord' 'Replay'}, % disable all buttons except Stop
        enable(A,0); enable(Q,0);
        enable(A('Stop'),1);
        GUIclosable(figh,0); % not okay to close GUI
        % color Play or PlayRec buttons
        if isequal('Play', Mode) || isequal('Replay', Mode),
            highlight(A('Play'),[0 0.7 0]);
        elseif isequal('PlayRecord', Mode),
            highlight(A('PlayRec'),[0.85 0 0]);
        end
end
% display GUI message, if any.
if ~isempty(gmess) && ~isequal(inf, gmess),
    GUImessage(figh,gmess,gmessmode);
end
% update the Experiment status info
EM = GUImessenger(figh, 'ExpInfo');
EXP = current(experiment);
reportstatus(EXP, EM);
figure(figh);
drawnow;

end

