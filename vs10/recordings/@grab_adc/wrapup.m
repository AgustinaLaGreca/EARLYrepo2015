function G = wrapup(G);
% grab_adc/wrapup - wrap up grab_adc object after finishing or stopping
%    G = wrapup(G) wraps up event grabbing by grab_adc object G.
%    The following steps are taken:
%       1. final call to grab_adc/oneshot to grab any pending events
%       2. store G's status @wrapup time in the G.Data.FinalStatus field.
%       3. delete the timer launched by grab_adc/start
%       4. change G's status to 'wrappedup'
%    Steps 3 & 4 are delegated to action/wrapup.
%
%    Note that no tranfer the data to G's receiving dataset is done here.
%    Data transfer is iniated by Dataset/getdata, which in turn calls upon 
%    @grab_adc/getdata to copy the data.
%
%   See also grab_adc/prepare, dataset/getdata, grab_adc/getdata.

eval(IamAt('indent')); % indent action method below
% 0. no download needed here, because all the calls do their own downloading and uploading.

% 1. final call to grab_adc/oneshot to grab any pending events
oneshot(G, 'sloppy'); % 'sloppy' = ignore status(G)

% 2. store G's status @wrapup time in the G.FinalStatus field.
download(G);
G = setdatabuf(G, 'FinalStatus', status(G));
upload(G);

% 3 & 4. delete the timer, change status, upload
wrapup(G.action); % download & upload are implicit





