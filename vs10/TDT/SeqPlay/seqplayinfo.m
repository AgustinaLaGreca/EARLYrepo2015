function y = SeqplayInfo(Check);
% SeqplayInfo - current status of Seqplay variables
%   SeqplayInfo returns a struct describing the current status of SeqPlay
%   variables.
%
%   SeqplayInfo(1) also checks whether the circuit is still loaded on the
%   TDT device.
%
%   See also SeqPlay, Seqplayinit, Sing.

if nargin<1, Check=0; end

y = private_seqPlayInfo; % info shared by seqplayXXX

if Check && ~isempty(y), % test if the circuit is really still there
    CI = sys3CircuitInfo;
    if isempty(CI), % circuit is not there anymore
        y = [];
    else, % check consistency of info
        [dum, circuitName] = fileparts(CI.CircuitFile);
        okay = isequal(CI.Device, y.Dev) ...
            && isequal(circuitName, y.Circuitname) ...
            && isequal(CI.Fsam, y.Fsam);
        if ~okay, y=[]; end
    end
end


