function E=resume(E);
% experiment/resume - resume experiment that is timed-out or finished.
%    E=resume(E) "touches" E, sets its state to 'Recording'  and makes it 
%    the current experiment.
%
%    E=resume(experiment()) resumes the Experiment that was most recently
%    the current experiment.
%
%    See also Experiment/iscurrent, Experiment/touch, Experiment/state,
%    Experiment/lastcurrentname.

if isvoid(E),
    E = find(E, lastcurrentname(E));
end

if iscurrent(E),,
    warning(['Experiment ''' name(E) ''' is already current. No need to resume.']);
elseif ~isvoid(current(experiment)),
    finish(current(experiment), '(resume)');
end
makecurrent(E);
E = touch(E);
E = state(E, 'Recording');
save(E, ['Resumed ' datestr(now)], 'new');
E = status(E,'Finished',[]);











