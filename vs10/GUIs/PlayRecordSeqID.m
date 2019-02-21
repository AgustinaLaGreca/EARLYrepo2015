function [ ds_info ] = PlayRecordSeqID(ds)
%play_record_seq_id After clicking on the play record button the user is
%asked to enter the cell number, extra stimulus type extension, pen depth
%and which electrode number
% ds contains cell and electrode information: Logical next one is retrieved from experiment status (status(Exp)), which holds
%last saved 
%   Detailed explanation goes here
figh= figure('color',[.75 .75 .75]);
set(figh,'position',[100 100 500 300]);
set(figh,'CloseRequestFcn', @local_cancel);

ui_struct.ds = ds;

txt = uicontrol('Style','text',...
        'Position',[50 209 75 20],...
        'String','Seq ID:','BackgroundColor',[.75 .75 .75]);
cell = uicontrol('Style','edit',...
        'Position',[150 210 50 20],...
        'String',num2str(ds.ID.iCell),'BackgroundColor','w');
ui_struct.cell = cell;

strup = char(hex2dec('2191'));
strdown = char(hex2dec('2193'));
cell_up = uicontrol('Style', 'pushbutton', 'FontName','Symbol','String',strup,...
    'Position', [165 234 20 20],...
    'Callback', @local_cell_up);    
cell_down = uicontrol('Style', 'pushbutton', 'FontName','Symbol','String',strdown,...
    'Position', [165 186 20 20],...
    'Callback', @local_cell_down);   
txt2 = uicontrol('Style','text',...
        'Position',[200 210 30 20],...
        'String','-','BackgroundColor',[.75 .75 .75]);
RecOfCell = uicontrol('Style','edit',...
    'Position',[230 210 50 20],...
    'String',num2str(ds.ID.iRecOfCell),'enable','off');
ui_struct.RecOfCell = RecOfCell;

txt3 = uicontrol('Style','text',...
        'Position',[280 210 30 20],...
        'String','-','BackgroundColor',[.75 .75 .75]);
StimType = uicontrol('Style','edit',...
    'Position',[310 210 50 20],...
    'String',ds.StimType,'BackgroundColor','w');
txt4 = uicontrol('Style','text',...
        'Position',[360 210 30 20],...
        'String','-','BackgroundColor',[.75 .75 .75]);
SubStimType = uicontrol('Style','edit',...
    'Position',[390 210 50 20],...
    'String','','BackgroundColor','w');
ui_struct.SubStimType = SubStimType;


St = status(current(experiment));

pen_depth_txt = uicontrol('Style','text',...
        'Position',[50 146 75 20],...
        'String','Pen Depth:','BackgroundColor',[.75 .75 .75]);
    
pen_depth = uicontrol('Style','edit',...
    'Position',[150 150 100 20],...
    'BackgroundColor','w','String',num2str(St.PenDepth));
micron_txt = uicontrol('Style','text',...
        'Position',[250 146 75 20],...
        'String','micron','BackgroundColor',[.75 .75 .75]);
ui_struct.pen_depth = pen_depth;

electrode_txt = uicontrol('Style','text',...
        'Position',[50 96 75 20],...
        'String','Electrode: ','BackgroundColor',[.75 .75 .75]);
% Electrode number zero is not allowed
if St.iPen == 0
    St.iPen = 1;
end
electrode = uicontrol('Style','edit',...
    'Position',[150 100 100 20],...
    'BackgroundColor','w','String',St.iPen);
ui_struct.electrode = electrode;

set(figh,'userdata',ui_struct);

OK = uicontrol('Style', 'pushbutton','String','OK',...
    'Position', [200 50 35 20],...
    'Callback', @local_ok);    
cancel = uicontrol('Style', 'pushbutton','String','Cancel',...
    'Position', [250 50 50 20],...
    'Callback', @local_cancel);    

uiwait(); % Then it is resumed after uiresume is executed
   ud = get(figh,'userdata');
   if isfield(ud,'return_struct')
       ds_info = ud.return_struct;
       delete(figh);
       return;
   elseif isfield(ud,'cancel')
       ds_info = [];
       delete(figh);
       return;
   end

end

function local_ok(source, event)
figh=gcf;
ud = get(figh,'userdata');
return_struct.iCell = str2num(get(ud.cell,'String'));
return_struct.iRecOfCell = str2num(get(ud.RecOfCell,'String'));
return_struct.StimType = ud.ds.StimType;
return_struct.SubStimType = get(ud.SubStimType,'String');
return_struct.pen_depth = str2num(get(ud.pen_depth,'String'));
return_struct.electrode = str2num(get(ud.electrode,'String'));
if isempty(return_struct.electrode)
    return_struct.electrode = 1;
end
ud.return_struct = return_struct;

set(figh,'userdata',ud);
uiresume();
end

function local_cancel(source, event)
figh=gcf;
ud = get(figh,'userdata');

ud.cancel = 1;

set(figh,'userdata',ud);
uiresume();

end

function local_cell_up(source,event)
    figh=gcf;
    
    ud = get(figh,'userdata');
    iCell = str2num(get(ud.cell,'String'));
    set(ud.cell,'String',num2str(iCell+1));
    if iCell+1 ~= ud.ds.ID.iCell
        set(ud.RecOfCell,'String','1');
    else
        set(ud.RecOfCell,'String',num2str(ud.ds.ID.iRecOfCell+1));
    end
end

function local_cell_down(source,event)
    figh=gcf;
    
    ud = get(figh,'userdata');
    iCell = str2num(get(ud.cell,'String'));
    if iCell-1 < 1
        warndlg('iCell must be larger then 0!','Error iCell');
        set(ud.cell,'String',num2str(iCell));
    else
        if iCell-1 ~= ud.ds.ID.iCell
            set(ud.RecOfCell,'String','1');
        else
            set(ud.RecOfCell,'String',num2str(ud.ds.ID.iRecOfCell+1));
        end
        set(ud.cell,'String',num2str(iCell-1));
    end
end
