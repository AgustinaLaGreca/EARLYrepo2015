function disp(Q)
% stimPresent/disp - DISP for stimPresent objects.
%
%   See stimpresent.

if numel(Q)>1,
    for ii=1:numel(Q),
        disp(Q(ii));
    end
    return
end
% single element from here

if isempty(Q.Fsam),
    disp('Void stimPresent object');
    return;
end

disp('&&&&&&&&&&&-stimPresentX-&&&&&&&&&&')
disp(struct(Q));
disp('&&&&&&&&&&&----&&&&-----&&&&&&&&&&')






