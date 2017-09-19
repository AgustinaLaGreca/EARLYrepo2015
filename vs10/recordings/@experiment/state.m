function s=state(E, s);
% experiment/state - get/set state of Experiment object
%   State(E) returns the state of Experiment object E.
%   The state is 'Void' for void experiments. Non-void experiments can be
%   either 'Initialized', 'Recording', or 'Finished'.
%
%   E = State(E, S) sets the state of E to S. The state cannot be set to 
%   'Initialized'; the other states are okay. Note that no implicit saving
%   is done.

if nargin<2, % get
    s = E.Status.State;
    if isempty(s),
        s = 'Void';
    end
else, % set
    if isvoid(E), 
        error('Cannot change the state of a void Experiment.');
    else,
        if ~ismember(s, {'Initialized' 'Recording' 'Finished'}),
            error(['Invalid experiment state ''' s '''.']);
        end
    end
    E.Status.State = s;
    s = E; % return arg
end









