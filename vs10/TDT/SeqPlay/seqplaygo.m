function y = SeqplayGo;
% SeqplayGo - trigger sequenced play over sys 3 device
%   SeqplayGo triggers the sequenced play previously specified
%   by SeqplayInit, SeqplayUpload and SeqplayList.
%
%   Type 'help Seqplay' to get an overview of sequenced playback.
%
%   See also Sing, SeqplayInit, SeqplayUpload SeqplayList, SeqplayHalt, Seqplaystatus.

SPinfo = private_seqPlayInfo; % info shared by seqplayXXX

% check SeqPlay status 
if isempty(SPinfo), error('Sequenced play not initialized. Call SeqPlayInit first.');
elseif ~isequal('listed', SPinfo.Status),
   error('No playlist specified for sequenced play. Call SeqplayList first.');
end

sys3setpar(-1, 'EndSample', SPinfo.Dev); % this aborts any ongoing D/A
sys3trig(1, SPinfo.Dev); % this resets the counters
sys3setpar(SPinfo.NsamTotPlay+1, 'EndSample', SPinfo.Dev); % this triggers conversion of correct # samples






