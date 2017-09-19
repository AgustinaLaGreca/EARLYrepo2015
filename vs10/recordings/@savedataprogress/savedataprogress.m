function S = savedataprogress(adDS, DT, dataviewer, figh, Name, You,DataviewerParam);
% savedataprogress - save dataset during data collection
%    S=savedataprogress(DS, DT, figh, 'XYZ', 'You') constructs an 
%    Action object S whose task it is to save a dataset for online data 
%    analysis by dataviewer during data collection. 
%    S is added to the list of Action of the figure
%    with handle figh under the name XYZ, and it will be finished as soon
%    as the Action object named You is finished. The dataset DS to be 
%    analyzed is assumed to have an address (in fact, it is okay to pass 
%    the address instead of the dataset). Every DT ms, an updated copy of 
%    DS is passed to a save function. 
%    At wrapup time, the data are once more updated and saved.
%
%    A special value for the Dataviewer is '-', in which case
%    savedataprogress calls Nope (=do nothing).
%
%    See also Action, dataset/getdata.

% if isequal('-', dataviewer), % no dataviewer specified; launch none
%     return;
% end

% persistent handle % use the same data analysis server as before...
% 
% try %...if the server is still alive
%     Execute(handle, 'pwd'); % trivial job
% catch % the RPC server is unavailable.
%     handle = [];
% end
% if isempty(handle)
%     ver = regexp(version,' ','split'); ver = ver{1}; dots = find(ver=='.'); ver = ver(1:dots(2)-1);
%     try
%         % Start Matlab Automation Server as Local Out-of-Process Server
%         handle = actxserver(['Matlab.Application.' ver]); 
%         Execute(handle, ['cd ' EarlyRootDir '\startupdir']);
%         Execute(handle, 'startup');
%     catch ME % data analysis server could not be launched
%         disp('Error while creating online analysis server:');
%         disp(ME.message);
%         for s=ME.stack(:).', disp(s); end
%         return;
%     end
% end 

if ~isa(adDS, 'address'), adDS = address(adDS); end
datasaver = cellify(fhandle('datasave')); % easier passing to feval

% DT = numel(dataviewer)*DT; % dynamic interval

Tlast = clock;
if ~isempty(DataviewerParam) 
    S = CollectInStruct(adDS, Tlast, DT, datasaver,dataviewer,DataviewerParam); % note that DT is not passed to Action constructor. Instead, the timer ..
else
    S = CollectInStruct(adDS, Tlast, DT, datasaver,dataviewer);
end
S = class(S,mfilename, action('initialized',[200 2000])); % ... interval is a fixed 200 ms; see ....
S = upload(S, figh, Name); % ...savedataprogress/oneshot why. The long start delay gives competing actions a chance 
afteryou(S, figh, You);