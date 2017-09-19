function Slope = JL_APthrSlope(Uidx, Slope);
% JL_APthrSlope - AP threshold slope in V/s
%    Slope = JL_APthrSlope(Uidx) returns the AP threshold of JL recording 
%    having unique recording index Uidx. The Slope value is retrieved from 
%    a database. The slope is used by JLgetrec for isolating APs.
%
%    JL_APthrSlope(Uidx, Slope) overrides the value in the database.
%
%    See also JLgetrec.

CFN = fullfile(processed_datadir, 'JL', 'JL_APthrSlope', mfilename);

if nargin<2, % get 
    
else, % set
end










