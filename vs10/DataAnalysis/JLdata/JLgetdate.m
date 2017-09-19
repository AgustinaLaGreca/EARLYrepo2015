function S = JLgetdate(S);
% JLgetdate - add abf file date info to JLbeatVar struct
%    S = JLgetdate(S) adds the following fields the the output S of
%    JLvarstats
%         abfdate: date of abf file as string, e.g., '17-Nov-2009 16:00:56'
%      abfdatenum: sortable date in numerical format, e.g. 7.3409e+005

N = numel(S);
abfdate = cell(1, N);
abfdatenum = cell(1, N);
for ii=1:N,
    fn = fullfile(S(ii).ABFdir, S(ii).ABFname);
    qq = dir(fn);
    abfdate{ii} = qq.date;
    abfdatenum{ii} = qq.datenum;
end
[S.abfdate] = deal(abfdate{:});
[S.abfdatenum] = deal(abfdatenum{:});

% reorder fieldnames of S
FNS = fieldnames(S);
Nf = numel(FNS);
S = orderfields(S, FNS([1:12 Nf-1 Nf 13:Nf-2]));









