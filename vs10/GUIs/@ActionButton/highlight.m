function highlight(A,clr)
% ActionButton/highlight - highlight/lowlight rendered actionButton
%    highlight(A,CLR) sets the action button's background color to CLR.
%    CLR is a 3-component RGB value, defaulting to [1 0 0], i.e. red.
%    
%    highlight(A,'default') resets the action button's background color to
%    its original value.
%
%    highlight(A,'random') sets the action button's background color to
%    a random RGB value.
%
%    For arrays A, each button A(k) is highlighted the same way.
%
%   See ActionButton, ActionButton/enable.

if nargin<2, clr = [1 0 0]; end

if numel(A)>1,
    for ii=1:numel(A),
        highlight(A(ii),clr);
    end
    return
end

if ischar(clr),
    [clr, Mess] = keywordMatch(clr, {'default' 'random'});
    error(Mess);
    switch clr,
        case 'default',
            clr = A.uicontrolProps.BackgroundColor;
        case 'random',
            clr = rand(1,3);
    end
end

set(A.uiHandles.Button, 'BackgroundColor', clr);









