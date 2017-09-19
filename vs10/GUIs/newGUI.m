function h = newGUI(Name, Title, CreatedBy, varargin);
% newGUI - open GUI figure
%    h = newGUI(GUIname, FigurebarTitle, CreatedBy, ...), opens a new 
%    figure to be used as a GUI and returns a handle to this figure.
%    GUIname is the name of the GUI, which is stored in the figure's
%    useradata (see getGUIdata). FigurebarTitle is the window name.
%    CreatedBy is a cell array, the first element of which is the name or
%    function handle of the mfile creating the GUI. Any remaining elements
%    of CreatedBy are the arguments passed to this creating mfile.
%    Thus a later call feval(CreatedBy{:}) should re-create the GUI.
%
%    See also GUIpanel, paramquery, gcg.

if ~ischar(Name) || ~isvarname(Name),
    error('Name input arg should be valid Matlab identifier (see ISVARNAME).');
end

if ~ischar(Title),
    error('Title input arg must be char string.');
end

if ~iscell(CreatedBy),
    error('CreatedBy input arg should be cell array.');
end

GSF = GUIsettings('GUIfigure');
% select those fields of GSF and varargin that match figure properties
DefaultProps = FullFieldnames(GSF,figureProperties, 'select');
ExplicitProps = FullFieldnames(struct(varargin{:}),figureProperties);
% merge default & explicit properties. The latter take precedence
figureProps = structJoin(DefaultProps, ExplicitProps);
% open the figure
h = figure('Name', Title, 'tag', Name, 'visible', 'off', figureProps);
h = h.Number;
setGUIdata(h,'GUIname', Name);
setGUIdata(h,'CreatedBy', CreatedBy);
setGUIdata(h,'OkayToClose', 1);



