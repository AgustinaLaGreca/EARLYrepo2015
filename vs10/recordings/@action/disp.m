function disp(A);
% action/DISP - DISP for action object
%    disp(A) displays content of A. Overload for more specific info.
%
%    See also collectdata, playaudio.

if isa(A,'dynamic') && hasbeenuploaded(A), A=download(A); end
disp(['Action object of type ''' class(A) '''; status=' A.Status '.']);

