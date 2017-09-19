function S = oneshot(S, force);
% savedataprogress/oneshot - elementary action of savedataprogress object
%    oneshot(S) gets an updated copy of the dataset of savedataprogress
%    object S and passes it to the datasave function of S.
%    This action is only performed if the last call to datasave is 
%    at least S.DT ms ago.
%
%    oneshot(S, 'force') ignores the interval since the last update and
%    invokes the datasaver and notification unconditionally.
%
%    See also savedataprogress.
persistent viewer_handle

if nargin<2, force=''; end
% The timer interval is a fast 200 ms, to prevent holding other processes;
% wrapping up can only start when *all* processes have finished.
% The "true" interval is S.DT ms, implemented manually below.

S = download(S); % need updated Tlast info
if isempty(force) && etime(clock, S.Tlast) < 1e-3*S.DT, return; end % skip if last true oneshot was too recent
S.Tlast = clock; 
% get dataset from address & update its data
DS = download(S.adDS); 
DS = getdata(DS);
EE = DS.ID.Experiment;
global datatest
datatest = DS;
% pass it to the dataviewer (including any trailing args). Make sure to
% prevent errors, as this would disrupt the handling of other, dependent, 
% action objects.

%     try
%         if isempty(viewer_handle) || ~ishandle(viewer_handle)
%             viewer_handle = figure;
%             set(viewer_handle,'Visible','On');
%         else
%             viewer_handle=clf(viewer_handle,'reset');
%             set(viewer_handle,'Visible','On');
%         end
%         dataview(DS,S.DataviewerParam{1},viewer_handle);
%         St = status(EE);
%         iRec = St.Ndataset+1;
%         iCell = max(1, St.iCell);
%         iRecOfCell = St.iRecOfCell+1;
%         stimType = DS.Stim.StimType;
%         fig_title = [EE.name ', REC: ' int2str(iRec) ', <' int2str(iCell) '-' ...
%             int2str(iRecOfCell) '-' stimType '>'];
%         ax=gca;
%         title(ax,fig_title,'fontsize', 12, 'fontweight', 'bold', 'interpreter', 'none');
%     end


try
    feval(S.datasaver{1}, DS);
    % notify the server of savedataprogress that an updated version of the
    % current dataset of experiment EE is available
%     Feval(S.handle,'datagraze', 0, name(EE));
    if ~isempty(force)
        pause(2);
        feval(S.datasaver{1}, DS, 'delete'); 
    end
    
catch ME,
    disp(['Error while handling savedataprogress object ''' name(S) ''' for online analysis:']);
    disp(ME.message);
    for s=ME.stack(:).', disp(s); end
end
upload(S);






