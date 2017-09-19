function disp(E);
% experiment/disp - disp for Experiment objects

if numel(E)>1,
    disp([sizeString(size(E)) ' Experiment object'])
elseif isvoid(E),
    disp('void Experiment object');
else,
    disp(['Experiment named ''' name(E) ''' of type ''' type(E) ''' started by experimenter ''' experimenter(E) ''', ' ...
        E.ID.Started '.']);
    if iscurrent(E),
        Cstr = ' *** Current Experiment *** ';
    else,
        Cstr = '';
    end
    disp(['   State: ' state(E) '.   ' Cstr]);
    disp(['   Definition last modified ', modified(E) '.']);
    disp(['       Status last modified ', statusmodified(E) '.']);
end









