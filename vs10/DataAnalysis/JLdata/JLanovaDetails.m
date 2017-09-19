function A = JLanovaDetails(Uidx);
% JLanovaStats - retrieve anova stats previously computed by JLanova2
%    W = JLanovaDetails(Uidx)
%    See JLanova2 -> local_cache for storage conventions.
%
%    See also JLanova2, JLanovaStats.

GBdir = fullfile(processed_datadir,  '\JL\NNTP');
load(fullfile(GBdir, ['Anova2_details_' num2str(Uidx)]), 'A');











