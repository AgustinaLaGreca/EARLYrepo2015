function A = Acellerator(A, MenuLabel, Key);
% ActionButton//Accelerator - attach accelerator key to action button
%    A = Acellerator(A, MenuLabel, Key) attaches an accelarator key to
%    actionbutton A. When A is rendered, the GUI will get a MenuItem having
%    the same callback as A, and accelerator key Key. The MenuItm will
%    placed at the uimenu labeled MenuLabel.
%
%    Accelerator only works with "simple" ActionButtons, i.e., non-toggling
%    ones.
%
%   See also ActionButton.

if nargout<1, error('Too few output args. Syntax is A=Acellerator(A,..)'); end

if ~ischar(MenuLabel), error('MenuLabel input argument must be char string.'); end
if ~isequal(1,numel(Key)) || ~isletter(Key),
    error('Key must be single letter.'); 
end
if iscellstr(A.String),
    error('Accelerators cannot be attached to toggling action buttons.');
end

A.Accel.MenuLabel = MenuLabel; 
A.Accel.Key = Key; 

A.Tooltip = ['[Ctr-' upper(Key) '] ' A.Tooltip];






