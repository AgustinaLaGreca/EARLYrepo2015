function F = save(E, LogStr, flag);
% experiment/save - save Experiment object.
%    save(E) saves the Experiment object E in filename(E) and extension 
%    '.ExpDef'. This is only allowed in the context of a running experiment, 
%    that is, when E is the current experiment and when the E.Computer is 
%    equal to Compuname().
%
%    save(E, LogStr) also saves char string Logstr to the E's logfile.
%
%    save(E, '', 'new') skips checking whether E is current. 
%
%    See also Datadir, experiment/filename, Experiment/addtolog.

if nargin<2,
    LogStr = '';
end
if nargin<3,
    flag = '';
end
checkfields(E); % check consistency of fieldnames

if isvoid(E), 
    error('Cannot save a void Experiment object.')
elseif ~isequal('new', flag) && ~iscurrent(E),
    error('E is not the current experiment and may therefore not be saved.');
elseif ~isequal('new', flag) &&  ~isequal(CompuName, computer(E)),
    error('Experiment object can only be saved on its native Computer.');
end

save([filename(E) '.ExpDef'],'E', '-mat');

if ~isempty(LogStr),
    addtolog(E, LogStr);
end











