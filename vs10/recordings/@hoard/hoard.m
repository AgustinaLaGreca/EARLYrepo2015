function H = hoard(Name, NsamChunk);
% hoard - construct hoard object for dynamic data storage
%    H = hoard creates a hoard object, which is an array that is accessed 
%    via an address (pointer) and that is safe from access by others.
%
%    See also eventdata, grabevents, dataset.

if nargin<1, Name = []; end
if nargin<2, NsamChunk = []; end

if isempty(Name),
    Name=num2str(RandomInt(1e9)); 
end


% properties of figure carrying the hoard
figProps = struct('visible', 'off', 'handlevisibility', 'off', ...
        'IntegerHandle', 'off', 'tag', Name, ...
        'userdata', 'dontclose'); % this saves it from popular aa and aaa figure killers
    
if ischar(Name), % create invisble, quasi indestructible, figure
    if ~isempty(findall(0, figProps)),
        error(['Hoard name ''' Name ''' is already in use.']);
    end
    figh = figure(figProps); 
    delete(allchild(figh)); % how cruel ..
    % data wil be stored as userdata of a uicontrol in figh:
    h = uicontrol(figh, 'visible', 'off', 'HandleVisibility', 'off', 'userdata',[], ...
        'tag', 'I_am_a_lonesome_hoard');   % leave a mark allowing H.figh to be recognized as hoard carrier
elseif ishandle(Name), % "reconstruct" existing hoard from its graphics handle
    figh = Name;
    h = allchild(figh);
end

Filename = fullfile(tempdir, ['hoard' Name]);
qq=memory; 
if isempty(NsamChunk), 
    NsamChunk = 1e6; % max # real doubles contained in one hoard
end
data = []; NchunkSaved = 0; NsamSaved = []; ScaleFactor = [];
set(h, 'userdata', CollectInStruct(data, NsamChunk, Filename, NchunkSaved, NsamSaved, ScaleFactor));
H = CollectInStruct(figh, h, Name);

H = class(H, mfilename);






