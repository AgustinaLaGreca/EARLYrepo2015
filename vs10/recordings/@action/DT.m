function DT=DT(A,DT);
% action/DT - set/get interval between background action shots
%    DT(A) returns the DT [ms] value that was specifed at contructor time.
%
%    A=DT(A,dt) sets interval to dt ms.
%
%    See also Action.

A = download(A,'sloppy');
if nargin<2, % get
    DT = A.DT;
else, % set
    A.DT = DT;
    if ~isempty(address(A)), DT = upload(A); end
end
    


