function click(h,dum,dumT,MouseButton); %#ok<INUSL>
% toggle/click - click callback of toggle
%    toggle(h,dum,T,'Left') is called by toggle button with handle h when
%    left-clicked. It sets string property of the button to the the next 
%    string from the list (rotating).
%
%    toggle(h,dum,T,'Right') is called by toggle button when right-clicked.
%    It enables/disables the state of the toggle button.

% get current state of toggle, i.e., the object in userdata clicked button
T = get(h,'userdata');
switch MouseButton,
    case 'Left', % regular click  : rotate to next string
        imatch = strmatch(T.Str0,T.StrArray, 'exact');
        if isempty(imatch), imatch = T.Nstr; end;
        inew = 1+rem(imatch,T.Nstr); % rotate
        T.Str0 = T.StrArray{inew};
        T=enable(T,1);
    case 'Right', % right-click" "enable/disable" button 
        T=enable(T,'swap');
end
% T has changed; re-render it
show(T);







