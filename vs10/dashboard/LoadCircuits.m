function LoadCircuits(CI);
% LoadCircuits - load RPvdS circuit according to info supplied by recordingInstructions
%    LoadCircuits(CI) loads RPvdS circuit(s) according to the second output 
%    argument of recordingInstructions.
%
%    See also recordingInstructions.

if ~isempty(CI.RX6_circuit),
    [Fsam, CycleUsage, Recycled] = sys3loadCircuit(CI.RX6_circuit, 'RX6', CI.RX6_Fsam/1e3,1); % last arg 1=recycle if possible
end

if ~isempty(CI.RP2_circuit),
    sys3loadCircuit(CI.RP2_circuit, 'RP2', CI.RP2_Fsam/1e3);  % last arg 1=recycle if possible
end




















