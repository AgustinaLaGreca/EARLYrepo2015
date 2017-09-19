function S = tk2xls(P, FN);
%  TK2XLS - convert TKpool info to XLS file
%    TK2XLS(P,'Foo'), where P is a TKpool struct array as returned by TKppol,
%    writes an overview of stimulus conditions to a file named Foo.xls.
%
%    See also TKpool, XLSWRITE.

FFN = fullfilename(FN, fullfile(EarlyRootDir,'sandbox'), '.xls');

S.ExpID = {};
S.Unit = {};
S.recID = {};
S.StimType = {};
S.SeqID = {};
S.ABF = {};
S.ISI = [];
S.BurstDur = [];
S.Nrep = [];
S.SPLi = []; S.SPLc = [];
S.Freqi = []; S.Freqc = [];
S.ITD = []; 
S.sane = [];
S.StimDetails = {};
for ii=1:numel(P),
    p = P(ii);
    Nrec = numel(p.Recs);
    try,
        DS = sgsrdataset(p.SGSRExpID, p.SGSRidentifier);
    catch,
        continue; % skip this one; it has no stimulus
    end
    [spl, fcar, dum, dum2] = SameSize(DS.SPL, DS.fcar, [1 1], DS.xval);
    Nchan = size(spl,2);
    Nmiss = Nrec-size(spl,1); % # missing conditions
    if Nmiss>0, % supply nans for missing conditions
        spl(end+(1:Nmiss),1:Nchan) = nan(Nmiss,Nchan);
        fcar(end+(1:Nmiss),1:Nchan) = nan(Nmiss,Nchan);
    end
    for irec=1:Nrec,
        S.ExpID = [S.ExpID; p.ExpID];
        S.Unit = [S.Unit; p.Pen];
        S.recID = [S.recID; p.recID];
        S.StimType = [S.StimType; DS.stimtype];
        S.SeqID = [S.SeqID; DS.SeqID];
        S.ABF = [S.ABF; p.Recs{irec}];
        S.ISI = [S.ISI; DS.RepDur];
        S.BurstDur = [S.BurstDur; unique(DS.Burstdur)];
        S.Nrep = [S.Nrep; DS.Nrep];
        S.SPLi = [S.SPLi; spl(irec,1)];
        S.SPLc = [S.SPLc; spl(irec,2)];
        S.Freqi = [S.Freqi; 1e-2*round(100*fcar(irec,1))];
        S.Freqc = [S.Freqc; 1e-2*round(100*fcar(irec,2))];
        if Nmiss>0 || ~isequal('NTD', DS.stimtype),
            itd = nan;
        else,
            itd = DS.xval(irec);
        end
        S.ITD = [S.ITD; round(itd)];
        S.sane = [S.sane; isequal(Nrec, DS.Nrec)];
        S.StimDetails = [S.StimDetails DS.info];
    end
end
local_xlswrite(S,FFN);


%======================================
function local_xlswrite(S,File)
if exist(File,'file'),
    delete(File);
end
ABC = char(double('A')+(0:25));
FNS = fieldnames(S).';
Ncol = numel(FNS);
xlswrite(File, FNS, ['A1:' ABC(Ncol) '1']);
for icol=1:Ncol,
    fn = FNS{icol};
    Nrow = numel(S.(fn));
    switch class(S(1).(fn)),
        case 'char'
            X = {S.(fn)}.'; % column cell array of strings
        otherwise % assume numerical ..
            X = cat(1,S.(fn)); % column array of doubles
    end
    xlswrite(File, X, [ABC(icol) '2:' ABC(icol) num2str(Nrow+1)]);
end














