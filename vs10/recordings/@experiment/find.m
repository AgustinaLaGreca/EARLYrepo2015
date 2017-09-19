function E = find(dum, ExpName);
% Experiment/find - find named experiment.
%    E = locate(experiment, 'RG09211') returns the experiment object E of
%    experiment RG09211. The first arg is a dummy, serving to make find a
%    Experiment method. A void experiment is returned when the named 
%    experiment cannot be found on this computer.
%
%    See Experiment/locate.

E = experiment(); % void exp
DD = locate(E, ExpName);
if ~isempty(DD) && exist(DD, 'dir'),
    FFN = fullfile(DD, [ExpName '.ExpDef']);
    warning('off', 'MATLAB:fieldConstructorMismatch');
    E = load(FFN, '-mat'); 
    warning('on', 'MATLAB:fieldConstructorMismatch');
    E = experiment(E.E);
    E = fixname(E); % fix misnamed experiment
    if isstruct(E), E = dum; end % don't load obsolete versions
end








