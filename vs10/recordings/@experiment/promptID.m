function E = promptID(E,figh);
% experiment/promptID - prompt for the ID of the current Dataset in this
% Experiment object's status
%   promptID(E,figh) - returns updated Experiment object and uses figh 
%   to extract the current Dataset's stimulus type from
%   promptID(E) - uses current GUI gcg to extract the current Dataset's
%   stimulus type form
%
%   See also Experiment/status, Experiment/save 

if nargin<2
    figh = gcg;
end

Pref = preferences(E);
if isequal(Pref.CheckID, 'No'), return; end; % no need for prompting

StimParam = getGUIdata(figh, 'StimParam');
stimType = StimParam.StimType;

StE = status(E); % status of experiment prior to this
iCell = max(1, StE.iCell);
iRecOfCell = StE.iRecOfCell+1;

% prompting
prompt={'dataset ID:'};
name='Confirm the ID for the current dataset.';
numlines=1;
defaultanswer={[num2str(iCell) '-' num2str(iRecOfCell) '-' upper(stimType)]};
okay=0;
while ~okay
    answer=inputdlg(prompt,name,numlines,defaultanswer);
    if isempty(answer) % user cancelled; use default values
        E=[];
        return; 
    end
    answer = regexp(answer{1},'-','split');
    if size(answer,2)<3
        waitfor(errordlg('Invalid input.','Input Error','replace'));
        continue;
    end
    newiCell = abs(round(str2double(answer{1}))); 
    newiRecOfCell = abs(round(str2double(answer{2})));
    newStimType = char(answer{3});
    if ~(issinglerealnumber(newiCell,'posinteger') && issinglerealnumber(newiRecOfCell,'posinteger'))
        waitfor(errordlg('Invalid Cell and/or Rec number.','Input Error','replace'));
        continue;
    end
    if ~isequal(newStimType,upper(stimType)) % this can't be changed
        waitfor(errordlg('Invalid Stim type.','Input Error','replace'));
        continue;
    end
    % the cell index may have been used before; 
    % find previous recordings from this cell
    iprev = find([StE.AllSaved.iCell]==newiCell); % indices of prev rec
    if isempty(iprev)
        newiRecOfCell = 1;
    else
        if ~(newiRecOfCell>max([StE.AllSaved(iprev).iRecOfCell]))
            waitfor(errordlg('Increase Rec number.','Input Error','replace'));
            continue;
        end
    end
    % passed all tests
    okay = 1;
end

% insertnote, if any
if size(answer,2)>3
    Note = ['% ' char(answer{4})];
    Head = ['% =========== note added ' datestr(now) '========='];
    Tail = repmat('=', size(Head)); Tail(1:2) = '% ';
    Note = strvcat(Head, Note, Tail);
    addtolog(E, Note);
end

E = status(E,'iCell',newiCell,'iRecOfCell',newiRecOfCell-1); % incremented in dataset.save
save(E);


