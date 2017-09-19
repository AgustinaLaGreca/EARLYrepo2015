function JLrefresh(icell_run);
% JLrefresh - refresh all JL caches
%    JLrefresh refreshes caches by calling all elementary analyses fro all
%    recordings.
%
%    JLrefresh(N) only refreshes the recordings of cell #N (icell_run index).
%
%    See also JLanova2, JLbeatvar.

if nargin<1,
    DB = JLdbase;
else,
    DB = JLdbase('icell_run', icell_run);
end

clear JLcycleCorrs

for irec=1:numel(DB),
    db = DB(irec)  % display to command line by way of progress report
    uix = db.UniqueRecordingIndex;
    JLbeatVar2(uix,0); % 0==no plot  - refreshes JLcycleCorrs output
    JLanova2(uix,0); % 0==no plot - refreshes output of JLanovaStats & JLanovaDetails
end











