function [Mess, Ignored, FN]=GUIfill(figh,S,doDisp, Except);
% GUIfill - restore stored string values and enable states of GUI paramqueries
%    Mess=GUIfill(h,S) sets the strings of edits and units of the paramqueries
%    of a GUI to the values previously collected by GUIgrab into struct S,
%    and restores the enable states of the queries.
%    The handle of the GUI is h. Any parameters contained in S that are not
%    present in the current GUI are ignored. Also, any illegal Unit value
%    (i.e. a value that is not in the current repertoire of the toggle) is
%    ignored. GUIfill(figh) returns a warning Mess that is only non-empty
%    when the name of the GUI does not match that of the GUIname stored in S.
%    If GUIfill is called with no output arguments, Mess is displayed as a
%    GUImessage.
%
%    Mess=GUIfill(h,S,1) forces displaying the message Mess, despite the
%    presence of an output argument.
%
%    [Mess, Ignored]=GUIfill(h,S) also returns a cellstr of Units that have
%    not been set due to incompatibility, e.g., because the stored string was
%    not a valid togglebutton string.
%
%    GUIfill(h,'Foo') retrieves the default values from file Foo in the GUI
%    defaults subdirectory used to store default GUI parameters of the
%    appropriate type. A GUI's "type" equals getGUIdata(h,'GUIname'); this
%    name is set by newGUI. Reports of file errors, etc, is done via the
%    output argument Mess if existent, displayed as a warning otherwise.
%
%    [Mess, Ignored, FN] = GUIfill(h,'?') prompts for the filename and 
%    then retrieves parameters from the file. The filename is returned in
%    output argument FN, which is empty if the user cancels.
%
%    GUIfill(h,'?',doDisp, {'Foo' 'Faa' ...}) retrieves all parameters 
%    Zxcept Foo, Faa, etc.
%
%    GUIfill(h,0) retrieves the parameter set of the GUI most recently
%    saved by GUIgrab. This typically means that this parameter set
%    has been been Checked & passed the test. 
%
%    GUIfill(h,0,doDisp, {'Foo' 'Faa' ...}) retrieves all parameters except Foo,
%    Faa, etc.
%
%    GUIfill(h,N, ...) retrieves the Nth parameter set counted backwards.
%
%    See also GUIval, GUIgrab, newGUI, GUImessage.

if nargin<3, doDisp=(nargout<1); end
if nargin<4, Except = {}; end

if ~istypedhandle(figh, 'figure'),
    error('First arg (figh) must be a handle to an existing figure.');
end
[Mess, Ignored, FN] = deal([]);
figure(figh);
GUIname = getGUIdata(figh,'GUIname');
if ischar(S), % a filename; retrieve struct S
    fn = S;
    DD = GUIdefaultsDir(GUIname, 'create'); % GUIdef dir of this GUI type; make sure it exists
    if isequal('?', fn),
        fn = uigetfile(fullfile(DD,'*.GUIdef'), ['Reading from file:  select ' GUIname ' parameter file']);
        if isequal(0,fn), FN=''; return; end % user canceled
        [~, fn, ~]= fileparts(fn); % strip off  extension
    end
    FN = FullFilename(fn, DD, '.GUIdef');
    if exist(FN,'file'),
        S = load(FN, '-mat'); S = S.grabStruct;
    else,
        Mess = (['GUIdefaults file ''' fn ''' not found.']);
        if doDisp,
            GUImessage(figh, Mess,'warn');
            FN = '';
        end
        return;
    end
elseif isnumeric(S), % most-recently used parameter set.
    M = S; % index of history set. M=0 is most recent one, M=1 is previous one, etc.
    DD = GUIdefaultsDir(GUIname, 'create'); % GUIdef dir of this GUI type; make sure it exists
    FN = FullFilename('most-recent', DD, '.GUIhistory');
    if exist(FN,'file'),
        S = load(FN,'-mat'); S = S.grabStruct;
        iS = length(S)-M;
        iS = max(1,iS);
        S = S{iS};
    else,
        return;
    end
end

if isempty(S), return; end % avoid needless trouble
% remove exceptions, if present
for ii=1:numel(Except),
    xp = Except{ii};
    if isfield(S,xp), 
        S = rmfield(S,xp);
    end
end

if isequal(GUIname, S.GUIname),
    Mess = '';
else,
    Mess = ['Filling GUI named ''' GUIname ''' using parameters obtained from GUI named ''' S.GUIname '''.'];
end
Ignored = {};

persistent GSQ; if isempty(GSQ), GSQ = GUIsettings('ParamQuery'); end

Q = getGUIdata(figh,'Query',[]); % array of queries.
% read queries one by one.
EnabStates = getFieldOrDefault(S, 'EnableStates_',[]); % enable states if present
for ii=1:numel(Q),
    q = Q(ii);
    EditStr = getFieldOrDefault(S,[q.Name 'Edit'],[]);
    he = getFieldOrDefault(q.uiHandles, 'Edit', []);
    UnitStr = getFieldOrDefault(S,[q.Name 'Unit'],[]);
    hu = getFieldOrDefault(q.uiHandles, 'Unit', []);
    enab = getFieldOrDefault(EnabStates, q.Name, []);
    if ~isempty(EditStr) && isUIcontrol(he), % impose string value on q 
        set(he,'string', EditStr, 'backgroundcolor', GSQ.EditOkayColor);
    end
    if isUIcontrol(he) && ~isempty(enab), % impose enable state on edit field
        set(he,'enable', enab);
    end
    if ~isempty(UnitStr) && iscellstr(q.Unit) && isUIcontrol(hu), % impose it on q
        T = get(hu,'userdata'); % toggle object
        [T, okay] = impose(T, UnitStr); % okay checks whether UnitStr is valid choice
        if okay,
            T = enable(T,1);
            show(T);
        else, Ignored = [Ignored, [q.Name 'Unit']];
        end
    end
    if iscellstr(q.Unit) && isUIcontrol(hu) && ~isempty(enab), % impose enable state on unit toggle
        set(hu,'enable', enab);
    end
end

if doDisp, % display message
    if ~isempty(Mess),
        GUImessage(figh, Mess,'warn');
    elseif exist('fn', 'var') && isempty(Ignored),
        GUImessage(figh, ['Parameters retrieved from file ''' fn '''.']);
    elseif exist('fn', 'var') && ~isempty(Ignored),
        GUImessage(figh, {['Parameters retrieved from file ''' fn ''''], ...
            'ignoring the following illegal Unit-values: ' Ignored{:}}, 'warn');
    end
end
