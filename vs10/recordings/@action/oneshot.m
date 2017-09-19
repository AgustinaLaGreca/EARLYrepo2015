function A=oneshot(A,varargin);
% action/oneshot - single action by action object
%    A=oneshot(A) performs the A's action a single time
%
%    Action is an abstract class, so action/oneshot must be overloaded 
%    for each concrete subclass. A call to the non-overloaded oneshot method
%    generates an error.
%
%    See also collectdata, playaudio.

% varargin{1}.Type
error(['action/oneshot is not overloaded for concrete action subclass ''' class(A) '''.']);



