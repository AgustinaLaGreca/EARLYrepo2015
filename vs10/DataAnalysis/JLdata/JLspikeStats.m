function PL = JLspikeStats(uix);
% JLspikeStats - retrieve spike stats previously computed by JLspikes
%    W = JLspikeStats
%    W = JLspikeStats(Uidx)
%    See JL_NNTP for storage conventions.
%
%    See also JLspikes, JLanovaDetails, JL_NNTP.

GBdir = fullfile(processed_datadir,  '\JL\NNTP');
load(fullfile(GBdir, 'PhaseLock.mat'), 'PL');
if nargin<1,
    PL = sortAccord(PL, [PL.UniqueRecordingIndex]);
else,
    PL = PL([PL.UniqueRecordingIndex]==uix);
end










