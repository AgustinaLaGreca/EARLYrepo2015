function G = wrapup(G);
% grab_comport/wrapup - wrap up grab_comport object after finishing or stopping
%    G = wrapup(G) wraps up event grabbing by grab_comport object G.
%    The following steps are taken:
%       1. store G's status @wrapup time in the G.Data.FinalStatus field.
%       2. perform G's CloseCall (this typically closes the serial connection)
%       3. delete the timer launched by grab_comport/start
%       4. change G's status to 'wrappedup'
%    Steps 3 & 4 are delegated to action/wrapup.
%
%    Note that no tranfer the data to G's receiving dataset is done here.
%    Data transfer is iniated by Dataset/getdata, which in turn calls upon 
%    @grab_comport/getdata to copy the data.
%
%   See also grab_comport/prepare, dataset/getdata, grab_comport/getdata.

eval(IamAt('indent')); % indent action method below
% 0. no download needed here, because all the calls do their own downloading and uploading.

% 1. store G's status @wrapup time in the G.FinalStatus field.
download(G);
G = setdatabuf(G, 'FinalStatus', status(G));
upload(G);

%2. call G's CloseCall, which typically closes the serial connection
feval(G.CloseCall{:});

% 3 & 4. delete the timer, change status, upload
wrapup(G.action); % download & upload are implicit





