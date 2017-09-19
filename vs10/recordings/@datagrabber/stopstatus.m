function G = stopstatus(G, X);
% datagrabber/stopstatus - get/set stopstatus of datagrabber object
%   G = stopstatus(G, X) sets G's stop  status to X.
%   St = stopstatus(G) returns G's stop  status.
%   datagrabber G. 
%
%   See also Grabevents/clear, action/clear.

if nargin<2, % get
    G = G.StopStatus;
else, % set
    G.StopStatus = X;
end

%   G = stopstatus(G, X) sets G's stop  statis to X.




