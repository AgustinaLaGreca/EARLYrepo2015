function T = JLcycleStats(Uidx);
% JLcycleStats - retrieve cycle stats previously computed by JLviewrec
%    T = JLcycleStats(Uidx)
%    See JL_NNTP for storage conventions.
%
%    See also JLwaveforms, JLviewrec, JL_NNTP.

persistent GGBB

GBdir = fullfile(processed_datadir,  '\JL\NNTP');

if isempty(GGBB),
    load(fullfile(GBdir, 'cycleStats'), 'T');
    GGBB = sortAccord(T, [T.UniqueRecordingIndex]);
end

if nargin<1,
    T = GGBB;
else,
    T = GGBB(ismember([GGBB.UniqueRecordingIndex], Uidx));
end








