function [L,N]=GUIactionList(figh, FN, flag);
% GUIactionList - list action objects uploaded to a GUI
%    GUIactionList(figh) returns a list addresses of action objects that
%    have been uploaded to the GUI with handle figh. An address is a struct
%    with fields get, put, rm as described in GUIaccess. The list is kept
%    and updated automatically each time when calling setGUIdata.
%
%    [L,N] = GUIactionList(figh) also returns the number of items N in the
%    list L.
%
%    GUIactionList(figh, 'Foo') adds the address of Foo to the action list of
%    this GUI. The address means GUIaccess(figh, 'Foo'). This type of call
%    should only be done by setGUIdata or else the bookeeping will be
%    screwed up.
%
%    GUIactionList(figh, 'Foo', 'remove') removes the address of Foo from
%    the action list of this GUI. This type of call should only be done by
%    setGUIdata or else the bookeeping will be screwed up.
%
%    GUIactionList(figh, 'init') creates an empty action list for figh, 
%    unless a list is already present, ion which case nothing is done. 
%    (only to be used by setGUIdata).
%
%    GUIactionList(figh, I, 'resort') changes the order of the current
%    action list stored in figh by replacing L by L(I). I must be a
%    permutation of 1:numel(L). This call affects the order of handling the
%    list members during a call to GUIaction.
%    (depricated)
%
%    See also GUIaction, setGUIdata, rmGUIdata.

if nargin<3, flag=''; end

L_em = emptystruct('FN','address'); % 0x0 struct with specified fields
L = getGUIdata(figh, 'GUIactionlist__', L_em);
if nargin>1,
    if isequal('init', FN), % initialize list, but only if none exists yet!
       if isequal(L,L_em),  % no list exists. Create one.
           setGUIdata(figh, 'GUIactionlist__', L_em); % nolist flag avoids endless GUIactions<->setGUIdata recursion
       end
       return;
    elseif isempty(L), ihit=[]; % not yet in list
    elseif isnumeric(FN),% wait
    else, ihit = strmatch(FN, {L.FN}, 'exact');
    end
    if isequal('remove', flag),
        L(ihit) = [];
        setGUIdata(figh, 'GUIactionlist__', L);
    elseif isequal('resort', flag),
        if ~isequal(sort(FN),1:numel(L)),
            error('In a resorting call to GUIaction, the index array I must be a permutation of 1:numel(L).');
        end
        setGUIdata(figh, 'GUIactionlist__', L(FN));
    elseif isempty(ihit), % add item
        item.FN = FN;
        item.address = GUIaccess(figh, FN);
        setGUIdata(figh, 'GUIactionlist__', [L item]);
    end
end
N = numel(L);







