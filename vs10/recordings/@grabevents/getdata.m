function D = getdata(G);
% grabevents/getdata - extract data data from Grabevents object
%   ED = getdata(G) returns Eventdata object ED holding the data events
%   collected by G sofar. This amounts to extracting the Data field of G 
%   with the pointer to the data, G.Data.EventTimes, replaced by the data 
%   themselves, get(G.Data.EventTimes). Before dereferencing the event
%   times, a copy of G is placed into Data.Supplier. As a final step, the 
%   resulting struct is cast into an Eventdata object. This prepares the 
%   data for transfer to the containing Dataset object, whether final or 
%   intermediate.
%
%   See also Grabevents, Dataset/getdata.

eval(IamAt);

Data = getdatabuf(G.datagrabber); 
Data.EventTimes = get(Data.EventTimes); % "dereference" the data

% cast into eventdata object
SupStatus = Data.FinalStatus;
if isempty(SupStatus), % G has not yet been visited by wrapup ..
    SupStatus = status(G); % we're in mid air, so get current status of G
end
D = eventdata(Data, SupStatus, G);







