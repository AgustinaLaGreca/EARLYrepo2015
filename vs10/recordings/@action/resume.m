function A=resume(A);
% action/resume - resume action of action object
%    resume(A) resumes the action of A that was previosuly stopped.
%
%    action/resume merely restarts the timer with handle thandle(A). If more
%    spefic work is needed for a certain subclass, the resume method
%    must be overloaded.
%
%    See also action/start, action/resume, collectdata, playaudio.

A=download(A);
if ~isequal('stopped', A.Status),
    error('Status of action object must be ''stopped'' in order to be resumed.');
end
start(A.thandle);
A=status(A,'resumed'); % implicit uploading here


