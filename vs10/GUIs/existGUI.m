function [eg,h] = existGUI(Name);
% existGUI - true when named GUI exists
%    existGUI('Foo') is true (1) if a GUI named Foo exists a a rendered
%    figure. Name must match exactly.
%
%    [E, H] = existGUI('Foo') also returns the graphics handle H of the 
%    named GUI if it exists, [] if it does not.
%
%    See also gcg, isGUI, newGUI.

% get all figure handles, also hidden ones
shh = get(0,'showhiddenhand'); % store to restore
set(0,'showhiddenhand', 'on');
hf = findobj(0,'type', 'figure');
set(0,'showhiddenhand', shh);

eg = 0; h = [];
for ii=1:numel(hf),
    hn(ii) = hf(ii).Number;
    eg = isGUI(hn(ii)) && isequal(Name,getGUIdata(hn(ii),'GUIname'));
    if eg,
        h = hn(ii);
        break;
    end
end


