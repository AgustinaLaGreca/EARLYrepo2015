function E=fixname(E);
% experiment/fixname - fix misnamed experiment
if isequal('RG11304', E.ID.Name) & isequal('Jolet', E.ID.Experimenter),
    E.ID.Name = 'RG11303'; 
end









