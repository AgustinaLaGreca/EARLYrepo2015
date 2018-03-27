function Q=draw(hp, Q, XY);
% ParamQuery/draw - draw paramquery object.
%   draw(hp, Q, X,Y) draws paramquery object Q in uipanel with handle hp at
%   position (X,Y) in pixels. This function is typically recursively called 
%   by GUIpiece/draw.
%
%   See ParamQuery, GUIpiece/draw, GUIpiece.
cpos=get(hp,'position');
X=XY(1); Y=cpos(4)-XY(2)-height(Q);
x=X; y=Y+0.25*height(Q);
% get widths & heights of individual parts of the query
[dum, Qwidth] = width(Q); 
[dum, Qheight] = height(Q);
drawEdit = ~isempty(Q.String); % do we need an Edit uicontrol?

% --prompt--
dy = 0.1*Qheight.Prompt; % correct vertical misplacement of Matlab text-style uicontrols
hprompt = uicontrol(hp, Q.uicontrolProps, 'position', [x y-dy Qwidth.Prompt Qheight.Prompt], ...
    'style', 'text', 'BackgroundColor', get(hp,'BackgroundColor'), ...
    'ForegroundColor', get(hp,'ForegroundColor'), ...
    'String', [Q.Prompt ' '], 'HorizontalAlignment', 'right', ...
    'TooltipString', Q.Tooltip);
x=x+Qwidth.Prompt;
% --edit--
if drawEdit,
    hedit = uicontrol(hp, Q.uicontrolProps, 'position', [x y Qwidth.Edit Qheight.Edit], ...
        'style', 'edit', 'BackgroundColor', get(hp,'BackgroundColor'), ...
        'ForegroundColor', get(hp,'ForegroundColor'), ...
        'String', Q.String, 'HorizontalAlignment', 'left');
    x=x+Qwidth.Edit;
end
% --unit--
un = Q.Unit;
if iscell(un), % toggle button
    if ~drawEdit, TT=Q.Tooltip; else, TT='Click toggle to select an option.'; end; % if Q is a toggle button w/o edit, button gets specific tooltip
    hunit = uicontrol(hp, Q.uicontrolProps, 'position', [x y Qwidth.Unit Qheight.Unit], ...
        'style', 'pushbutton','BackgroundColor', get(hp,'BackgroundColor'), ...
        'ForegroundColor',get(hp,'ForegroundColor') , ...
        'String', un{1},'HorizontalAlignment', 'right', 'Tooltip', TT);
      hunit = double(hunit);
      betoggle(hunit,un); % make it a toggling button
 
else, % fixed string
    dy = 0.1*Qheight.Unit; % correct vertical misplacement of Matlab text-style uicontrols
    hunit = uicontrol(hp, Q.uicontrolProps, 'position', [x y-dy Qwidth.Unit Qheight.Unit], ...
        'style', 'text', 'BackgroundColor', get(hp,'BackgroundColor'), ...
        'ForegroundColor', get(hp,'ForegroundColor'), ...
        'String', un, 'HorizontalAlignment', 'right');
end


% add Query info (including handles!) to userdata of GUI figure
figh = parentfigh(hp);
figh = figh.Number;
% Q.uiHandles.Unit = hunit;
Q.uiHandles.Unit = double(hunit);
Q.Parent = hp;
if drawEdit, 
    Q.uiHandles.Edit = double(hedit); 
    hct  = findobj(figh, 'tag', 'QueryContextMenu');
    set(hedit, 'userdata', Q, 'UIContextMenu', hct, ...
        'ButtonDownFcn', {@PassYourselfToContextMenu, ParamQuery});
end
Query = getGUIdata(figh, 'Query', []);
Query = [Query Q];
setGUIdata(figh, 'Query', Query);


%               Name: 'StartFrequency'
%             Prompt: 'Start:'
%             String: '12000 12001'
%               Unit: 'Hz'
%         Constraint: 'rreal/positive'
%            Tooltip: 'First frequency of series.'
%             MaxNum: 2
%     uicontrolProps: [1x1 struct]
% % % % % %         Units: 'pixels'
% % % % % %      FontSize: 8
% % % % % %      FontName: 'MS Sans Serif'
% % % % % %     FontUnits: 'points'

% %% lineCallback for toggling button
%     function lineCallback(src,~,arg)
%         Str = src.String;
%         n = numel(arg);
%         for i=1:n
%             if (strcmp(Str,arg{i})) && (i<n)
%                 src.String = arg{i+1};
%             elseif (strcmp(Str,arg{i})) && (i == n)
%                 src.String = arg{1};
%             end;
%         end;
%         
% 




