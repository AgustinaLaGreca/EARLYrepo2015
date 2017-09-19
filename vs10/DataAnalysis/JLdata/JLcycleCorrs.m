function T = JLcycleCorrs(Uidx);
% JLcycleCorrs - retrieve cycle correlation previously computed by JLbeatvar2
%    T = JLcycleCorrs(Uidx)
%    See JLbeatvar2 for storage conventions.
%
%    See also JLcycleStats, JLwaveforms, JLviewrec, JL_NNTP.

persistent GGBB

GBdir = fullfile(processed_datadir,  '\JL\NNTP');

if isempty(GGBB),
    load(fullfile(GBdir, 'cycleCorrs'), 'CC');
    GGBB = sortAccord(CC, [CC.UniqueRecordingIndex]);
end

if nargin<1,
    T = GGBB;
else,
    T = GGBB(ismember([GGBB.UniqueRecordingIndex], Uidx));
end








