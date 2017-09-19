function y = SeqplayHalt;
% SeqPlayHalt - immediately halt sequenced play D/A
%   SeqplayHalt triggers the seqplay circuit on the Sys3 device to stop D/A conversion.
%   The circuit is immediately halted and reset.
%
%   See also Seqplay, SeqPlayGo.

SPinfo = private_seqPlayInfo;
if ~isstruct(SPinfo),
    warning('No seqplay initialized');
    return;
end
sys3setpar(-1, 'EndSample', SPinfo.Dev); % this aborts any ongoing D/A
sys3trig(1,SPinfo.Dev); % this resets the counters
sys3halt(SPinfo.Dev);


