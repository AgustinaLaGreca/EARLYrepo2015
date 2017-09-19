function n=computer(E);
% experiment/computer - computer of Experiment object
%   Computer(E) returns the computer of Experiment object E.

n = E.ID.Computer;
if isempty(n),
    n = '';
end






