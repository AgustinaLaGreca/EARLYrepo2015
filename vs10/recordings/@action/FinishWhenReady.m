function A=FinishWhenReady(A);
% action/FinishWhenReady - test if action is finished; if so, stop timer
%    FinishWhenReady(A) tests whether the action spawned by A is finished. 
%    If so, the A's timer is stopped, and A's status is updated to
%    'finished'. Note that this regular ending of A's action is different
%    from the externally triggered ending by action/stop.
%
%    See also action/isready action/start, action/stop.
%    

A=download(A);
if isready(A),
    stop(A.thandle);
    A=status(A,'finished'); % implicit uploading here
end


