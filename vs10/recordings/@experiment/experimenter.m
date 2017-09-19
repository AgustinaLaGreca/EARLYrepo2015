function n=experimenter(E);
% experiment/experimenter - Experimenter of experiment object
%   Experimenter(E) returns the experimenter of Experiment object E.
n = E.ID.Experimenter;
if isempty(n), 
    n = ''; 
else,
    n(1) = upper(n(1)); % it's a proper name
end







