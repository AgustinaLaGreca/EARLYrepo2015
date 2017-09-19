function N = Nevent(ED);
% eventdata/Nevent - number of events stored in eventdata object
%   Nevent(ED) returns the number of events stored in eventdata object ED.
%
%   See also eventdata/eventtimes.

N = numel(ED.Data.EventTimes);







