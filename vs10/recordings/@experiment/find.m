function E = find(dum, ExpName);
% Experiment/find - find named experiment.
%    E = locate(experiment, 'RG09211') returns the experiment object E of
%    experiment RG09211. The first arg is a dummy, serving to make find a
%    Experiment method. A void experiment is returned when the named 
%    experiment cannot be found on this computer.
%
%    See Experiment/locate.

% Marta Oct18
% Introduced to change the Resume Experiment input GUI from typing name
% experiment to selecting from files explorer
% To revert this addition, delete until next commented out section and
% uncomment the aforementioned section
%% Section to be removed/commented out
E = experiment(); % void exp
[D, name] = fileparts(ExpName); % remove XX
if ~isempty(D)
    DD = fullfile(D,name);
    ExpName = name;
else
    DD = locate(E, ExpName);
end
%%  Uncomment this section for previous input gui
% E = experiment(); % void exp
% DD = locate(E, ExpName);
%%
if ~isempty(DD) && exist(DD, 'dir'),
    FFN = fullfile(DD, [ExpName '.ExpDef']);
    warning('off', 'MATLAB:fieldConstructorMismatch');
    try
        E = load(FFN, '-mat');
        warning('on', 'MATLAB:fieldConstructorMismatch');
        E = experiment(E.E);
        E = fixname(E); % fix misnamed experiment
        if isstruct(E), E = dum; end % don't load obsolete versions
    catch
        msg = ['Folder ', ExpName, ' is not a valid experiment folder.'];
        title = 'File Error';
        errorh = errordlg(msg,title);
        uiwait(errorh);
    end
end








