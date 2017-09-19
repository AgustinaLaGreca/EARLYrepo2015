function S = ZWOAEsummary(FN, idata);
% ZWOAEsummary - extract rel. pars for summary of ZWOAE rec
%   ZWOAEsummary('foo') summarizes pars in file foo.
%
%   ZWOAEsummary('foo', 15) summarizes pars in file foo015.
%
%
%  See also getZWOAEdata.

if length(idata)>1, % multiple sets: recursive handling
    S = [];
    for ii=1:length(idata),
        try,
        S = [S, ZWOAEsummary(FN, idata(ii))];
        catch, warning(['Nonexistent data set #'  num2str(idata(ii)) '''.']);
        end
    end
    return;
end

%-----single idata from here-----------

% read data
D = getZWOAEdata(FN, idata, '-nosig');
% extract important stimulus params
sp = D.stimpars;
iF1 = [sp.N1~=0] & [sp.L1>0] ;
iF2 = [sp.N2~=0] & [sp.L2>0] ;
F1 = sp.F1(iF1).'; % lower freq(s) in Hz
F2 = sp.F2(iF2).'; % upper freq(s) in Hz

N1 = sum(sp.N1);
N2 = sum(sp.N2);
f1 = mean(F1);
f2 = mean(F2);
df1 = mean(diff(F1)); if N1<2 df1 = 0; end
df2 = mean(diff(F2)); if N2<2 df2 = 0; end
df = max(df1,df2);
idata = D.index;

S = CollectInStruct(idata, N1, N2, f1, f2, df1, df2, df);



