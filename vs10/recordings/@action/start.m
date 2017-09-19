function A=start(A);
% action/start - start action of action object
%    A=start(A) starts the action that A does for a living. 
%    The timer with handle A.thandle is created (if none is present yet), 
%    and started (see action/createtimer & action/prepare), and A's status 
%    is set to 'started'.
%
%    For more complex type of actions, start may be overloaded for the
%    subclass at hand. Make sure to create the timer!
%
%    See also action/prepare, collectdata, playaudio.

eval(IamAt);
A = download(A);
if isempty(A.thandle), % attach timer
    A.thandle = createtimer(A);
end
A.Status = 'started'; 
upload(A);
start(thandle(A)); 



