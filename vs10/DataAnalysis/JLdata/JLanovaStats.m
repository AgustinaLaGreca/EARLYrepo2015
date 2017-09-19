function Astats = JLanovaStats(uix);
% JLanovaStats - retrieve anova stats previously computed by JLanova2
%    W = JLanovaStats
%    W = JLanovaStats(Uidx)
%    See JL_NNTP for storage conventions.
%
%    See also JLanova2, JLanovaDetails, JL_NNTP.

GBdir = fullfile(processed_datadir,  '\JL\NNTP');
load(fullfile(GBdir, 'Anova2_Stats'), 'Astats');
if nargin<1,
    Astats = sortAccord(Astats, [Astats.UniqueRecordingIndex]);
else,
    Astats = Astats([Astats.UniqueRecordingIndex]==uix);
end










