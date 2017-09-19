function D = getdata(G);
% grab_comport/getdata - extract data data from grab_comport object
%   ED = getdata(G) returns Eventdata object ED holding the data events
%   collected by G sofar. This amounts to retrieving the DataBuf of G and
%   replacing any datapointers (hoard objects) within DataBuf by real data.
%   This prepares the data for transfer to the containing Dataset object, 
%   whether final or intermediate.
%
%   See also grab_comport, Dataset/getdata.

eval(IamAt);

Data = getdatabuf(G.datagrabber); 
% "dereference" the data
Data.Samples = get(Data.Samples); 
Data.Times = get(Data.Times); 

% cast into eventdata object
SupStatus = Data.FinalStatus;
if isempty(SupStatus), % G has not yet been visited by wrapup ..
    SupStatus = status(G); % we're in mid air, so get current status of G
end
D = comport_data(Data, SupStatus, G);






