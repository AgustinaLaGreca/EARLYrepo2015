function T=enable(T,State); 
% toggle/enable - enable/disable toggle object
%    T=enable(T) or T=enable(T,1) enables T. 
%
%    T=enable(T,0) disables the togglebutton.
%
%    T=enable(T,'swap') swaps the togglebutton.
%
%    Enabling/disabling is triggered by right-clicking the button. Note 
%    that enable does not render the change - it only changes T.
%
%    See also betoggle, toggle/click, toggle/enableState.

if nargin<1, State = 1; end
if nargout<1, error('Too few output args.'); end

if isequal('swap', lower(State)),
    State = 1-T.EnableState;
end

T.EnableState = State;





