function E = openExpDef(FN);
%   openExpDef - display experiment definition from file.
%      openExpDef('Foo.ExpDef') displays the definition of Experiment Foo
%      as stored in file Foo.ExpDef.
%      
%   See also Experiment/load, Experiment/edit.

[Path, Nam] = fileparts(FN);
E = find(experiment, Nam);
disp(E)





