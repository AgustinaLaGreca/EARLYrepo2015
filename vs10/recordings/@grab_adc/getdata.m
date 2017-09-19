function D = getdata(G);
% grab_adc/getdata - extract data data from grab_adc object
%   ED = getdata(G) returns Eventdata object ED holding the data events
%   collected by G sofar. This amounts to extracting the Data field of G
%   with the Data.Samples subfield, which contains a pointer to the data,
%   and replacing it by the data themselves, get(G.Data.Samples). Before 
%   dereferencing the samples, a copy of G is placed into Data.Supplier. 
%   As a final step, the resulting struct is cast into an adc_data object. 
%   This prepares the data for transfer to the containing Dataset object, 
%   whether final or intermediate.
%
%   See also grab_adc, Dataset/getdata.

eval(IamAt);
Data = getdatabuf(G.datagrabber); % extract the data 
Data.Samples = getuserdata(Data.Samples); % dereference the data

% cast into adc_data object
SupStatus = Data.FinalStatus;
if isempty(SupStatus), % G has not yet been visited by wrapup ..
    SupStatus = status(G); % we're in mid air, so get current status of G
end
D = adc_data(Data, SupStatus, G);







