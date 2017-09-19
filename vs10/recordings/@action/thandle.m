function h=thandle(A,h);
% action/thandle - get/set thandle of action object
%    thandle(A) returns the current thandle of action object A. If A is 
%    dynamic, it is downloaded before returning its thandle.
%
%    A=thandle(A,h) sets the thandle of A to h. If A is dynamic, it is
%    downloaded before, and uploaded after its thandle has been set.
%
%    See also collectdata, playaudio.

HBU=hasbeenuploaded(A);
if HBU, A = download(A); end
if nargin<2, % get
    h = A.thandle;
else, % set
    A.thandle = h;
    if HBU, upload(A); end
    h = A;
end


