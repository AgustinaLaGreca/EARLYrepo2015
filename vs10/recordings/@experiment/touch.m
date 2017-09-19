function E = touch(E);
% experiment/touch - update statusmodified date of experiment.
%    touch(E) set the StatusModifed property of E to datestr(now). This
%    serves to overcome the restriction that Status, which is the regular
%    way of changing E's Status properties, may only be done when E is
%    current. E is saved after changing the property.
%
%    See also Experiment/iscurrent, Experiment/resume.

Dstr = datestr(now);
E.Status.StatusModified = Dstr;
addtostatuslog(E, ['touched ' Dstr]);
save(E,'', 'new');









