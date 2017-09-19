function S = oneshot(S, force);
% showdataprogress/oneshot - elementary action of Showdataprogress object
%    oneshot(S) gets an updated copy of the dataset of Showdataprogress
%    object S and passes it to the dataviewer function of S. This action is
%    only performed if the last call to the dataviewer is at least S.DT ms
%    ago.
%
%    oneshot(S, 'force') ignores the interval since the last update and
%    invokes the dataviewer unconditionally.
%
%    See also Showdataprogress.

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
% pass it to the dataviewer (including any trailing args). Make sure to
% prevent errors, as this would disrupt the handling of other, dependent, 
% action objects.
try,
    feval(S.dataviewer{1}, DS, S.dataviewer{2:end});
catch ME,
    disp(['Error while handling showdataprogress object ''' name(S) ''':']);
    disp(ME.message);
    for s=ME.stack(:).', disp(s); end
end
upload(S);






