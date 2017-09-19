function bc=isbeingcollected(D);
% Dataset/isbeingcollected - true dataset objects that are being collected
%   isbeingcollected(D) returns TRUE if D is a dataset whose data are still
%   being "grabbed", FALSE otherwise. 
%
%   See Dataset/getdata.

bc = ~isequal('complete', D.Status) && ~isequal('interrupted', D.Status);


