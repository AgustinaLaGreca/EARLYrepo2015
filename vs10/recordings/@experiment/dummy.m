function E=dummy(~, name);
% experiment/dummy - named dummy experiment 
%    dummy(experiment(), 'Foo') returns a dummy experiment named Foo.
%    For an application, see adc_data/folder.

E = experiment(); % make sure we return a dummy
E.ID.Name = name;
E.ID.Experimenter = 'cheater'; % this is to make isvoid(E) false.







