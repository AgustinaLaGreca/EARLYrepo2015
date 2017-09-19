function [Y,N] = recGUIpanel(E, kw, varargin);
% recGUIpanel - GUIpanel for recording settings
%   Qnames = recGUIpanel(Exp) returns GUIpanel to be rendered on
%   dashboard GUI. Exp is the current Experiment. P contains a number of 
%   toggles for enabling/disabling individual recording channels as
%   specified during the Experiment definition. Qnames is a cell array
%   containing the names of the toggles (which are ParamQuery objects).
%
%   See also StimGUI, Experiment.

persistent RecGUInames
if nargin<2, kw='stimGUIpanel'; end
switch kw,
    case 'stimGUIpanel', % recGUIpanel(E, 'stimGUIpanel')
        [Y,N] = local_panel(E, varargin{:});
        RecGUInames = N;
    case 'savesettings' % recGUIpanel(E, 'savesettings', figh)
        figh = varargin{1};
        GG = GUIgrab(figh);
        S = structpart(GG, local_appendUnit(RecGUInames));
        S.GUIname = GG.GUIname;
        putcache(mfilename,10, E.ID.Started,S);
    case 'restoresettings', % recGUIpanel(E, 'restoresettings', figh)
        figh = varargin{1};
        S = getcache(mfilename, E.ID.Started);
        if isempty(S),
            S = structpart(E.Recording.Source, RecGUInames);
            S = local_appendUnit(S);
            if isempty(S), clear S; end
            S.GUIname = getGUIdata(figh, 'GUIname');
        end
        if ~isempty(S), 
            try, GUIfill(figh, S); % not worth crashing
            end
        end;
    otherwise,
        error(['Invalid keyword ''' kw '''.']);
end

%================================================================
%================================================================
function [RP, FNS] = local_panel(E, varargin);
XArgs = {'' 'click to toggle recording'};
RP = GUIpanel('RecPanel', 'recording', varargin{:});
Src = E.Recording.Source;
if isempty(Src),
    FNS = {};
    RP = marginalize(RP, [100 60]);
    return;
end
% Each recording source of the experiment gets its own paramquery
FNS = fieldnames(Src);
Prompt_texts = cellstr(strjust(strvcat(FNS))); % right-justified version of prompts
for ii=1:numel(FNS),
    fn = FNS{ii};
    prmpt = [Prompt_texts{ii} ':'];
    PQ = ParamQuery(fn, prmpt, '', local_strip({Src.(fn).DataType}), XArgs{:});
    if ii==1, pos={'below' [7 0]}; else, pos = {'alignedwith' [0 -5]}; end
    RP = add(RP,PQ, pos{:});
end

function C=local_strip(C);
if isequal([],C{1}) || ~isempty(strmatch('-',lower(C), 'exact')),
    C = {local_dashes};
end
C = fliplr(unique([C local_dashes]));

function D=local_dashes;
D = '     -     ';

function C=local_appendUnit(C);
% append 'Unit' to all cells of C or fieldnames of struct C
if isstruct(C),
    S = C;
    FN = fieldnames(S);
    for ii=1:numel(FN),
        fn = FN{ii};
        C.([fn 'Unit']) = S.(fn);
    end
elseif isempty(C),
    C = {};
else,
    for ii=1:numel(C),
        C{ii} = [C{ii} 'Unit'];
    end
end


