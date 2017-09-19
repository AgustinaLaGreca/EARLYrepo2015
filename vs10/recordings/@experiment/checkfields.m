function checkfields(E);
% experiment/checkfields - Test data structure of Experiment object
%   checkfields(E) tests if Experiment E has the right data format, by
%   comparing it to a void Experiment object. An error is thrown when there
%   is an inconsistency.
%
%   See also Experiment.

E0 = experiment(); % reference: void object
FNS = fieldnames(E0);
for ii=1:numel(FNS),
    fn = FNS{ii};
    ref = E0.(fn);
    tst = E.(fn);
    if isequal('Version', fn),
        if ~isequal(ref,tst),
            error('Obsolete Version number pf Experiment object.');
        end
    elseif isstruct(ref),
        if ~isequal(fieldnames(ref),fieldnames(tst)),
            disp('===ref=====')
            disp(fieldnames(ref).');
            disp('===tst=====')
            disp(fieldnames(tst).');
            warning(['Incorrect data structure of field ''' fn ''' of Experiment object.']);
        end
    end
end










