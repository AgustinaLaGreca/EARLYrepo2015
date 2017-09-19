function clear(G);
% grab_comport/clear - clear grab_comport object and its data from memory
%   cleardata(G) removes Grabevent object G from history by offloading it
%   and by clearing from memory the data previously grabbed by G.
%   The offloading part is delegated to Action/clear.
%
%   See also action/clear, grab_comport/wrapup, dataset/getdata.

eval(IamAt('indent')); % indent action method below

% clear data from memory (remember that EventTimes are in a hoard object)
clear(getdatabuf(G, 'Samples'));
clear(getdatabuf(G, 'Times'));

% checking & offloading 
clear(G.action);





