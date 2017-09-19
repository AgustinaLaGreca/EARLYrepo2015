% Ruedata processing collection 
[Z dt Dur]=readrueData('A',1);
Z = Z/std(Z);
SigA=MapZ2datadistr(Z);
f1; hA = cdfplot(Z); hold on

[Z dt Dur]=readrueData('B',1);
Z = Z/std(Z);
SigB=MapZ2datadistr(Z);
f1; hB = cdfplot(Z); hold on

[Z dt Dur]=readrueData('Z',1);
Z = Z/std(Z);
SigZ=MapZ2datadistr(Z);
f1; hZ = cdfplot(Z); hold on

f1; hN = cdfplot(randn(1,1e4)); hold on

xlim([-5 5]);
set(hA, lico(1))
set(hB, lico(2))
set(hZ, lico(3))
set(hN, 'color', 'k', 'linewidth', 2);
legend({'A' 'B' 'Z' 'Normal'}, 'location','northwest')
set(gcf,'units', 'normalized', 'position', [0.28 0.577 0.665 0.322])

% Za = SigA(rand(1,numel(Z)));
% h = cdfplot(Z);
% set(h,'color', 'r', 'linestyle',':');

%==============================================
%   autocorr & ruecorr
RueTriplot('p90206Ac', 'p90206Bc'); RueTriplot('p90604Ac', 'p90604Bc'); RueTriplot('p90305Ac', 'p90305Bc'); RueTriplot('p90220Ac', 'p90220Bc')
RueTriplot('p90619Ac', 'p90619Bc');
RL = Ruelist
for ii=1:round(numel(RL)/2), for jj=ii+1:20, RL{ii}, RL{jj}, RueCrossShuf(RL{ii}, RL{jj}); end ; end

%=======population stuff=============================
RL = Ruelist; % all recordings
RLc = RL(~cellfun(@isempty, regexp(RL, 'c'))) % only the contra recs
S = RueCorrPop(RLc); % comparing across all pairs, identicals included
save D:\Data\RueData\IBW\Including_Zero_dB\FCCA\whole_audiogram_xcorr.mat S
load D:\Data\RueData\IBW\Including_Zero_dB\FCCA\whole_audiogram_xcorr.mat 
%-same analysis on interleaved c+i recordings
RL = Ruelist; % all recordings
RLc = RL(~cellfun(@isempty, regexp(Ruelist, 'c'))); % only the contra recs
RL2 = cellfun(@(str)str(1:length(str)-1), RLc, 'UniformOutput', false) % remove trailing 'c'
S2 = RueCorrPop(RL2); % comparing across all pairs, identicals included
save D:\Data\RueData\IBW\Including_Zero_dB\FCCA\whole_audiogram_xcorr.mat S2 -append
load D:\Data\RueData\IBW\Including_Zero_dB\FCCA\whole_audiogram_xcorr.mat 
%
qnonpair = ~[S.isSame] & ~[S.isPair]; qpair = [S.isPair]; qSame = [S.isSame];
% ===check effect of excluding the diagonal terms for pairs & non-pairs
plot([S(qnonpair).rho_nd], [S(qnonpair).rho_diag], '.');
xplot([S(qpair).rho_nd], [S(qpair).rho_diag], 'o');
axis equal; xlim(xlim);ylim(ylim); xplot([-1 1], [-1 1], 'k');
xlabel('non-diagonal xcorr', 'fontsize', 12); ylabel('diagonal xcorr', 'fontsize', 12);
legend({'Non pairs' 'Pairs'}, 'location' ,'northwest')
% ===distribution of non-diagonal autocorr with cells
NSame = sum(qSame); SameName = {S(qSame).FN1}; iSame = find(qSame);
%stem(1:NSame, [S(qSame).rho_nd]);
Ss = S(qSame); rn_ = [Ss.rho_nd]; st_ = [Ss.std_nd]; tt_ = logical([Ss.ttest_nd]);
errorbar(1:NSame, rn_, st_, 'linestyle', 'none', 'marker', 's', 'markerfacecolor', 'w');
xplot(find(tt_), rn_(tt_), 'linestyle', 'none', 'marker', 's', 'markerfacecolor', 'b')
ylim([0 1]); 
for ii=1:NSame, text(ii,S(iSame(ii)).rho_nd+0.05, SameName{ii}, 'rotation', 75); end
set(gca,'xtick',[]); title('Non-diagonal autocorrelation (whole audiogram)');
DeafOnes  = {Ss(~tt_).FN1}; Ndeaf = numel(DeafOnes);
figure; for ii=1:Ndeaf, subplot(3,2, ii); RueceptiveField(DeafOnes{ii}); end
% === distr of pairs vs nonpairs
Sp = S(qpair); Snp = S(qnonpair); Npair = numel(Sp); Nnonpair = numel(Snp);
Sp_dis = Sp([Sp.tt_both_1]); Snp_dis = Snp([Snp.tt_both_1]); Npair_dis = numel(Sp_dis); Nnonpair_dis = numel(Snp_dis);
rrr = linspace(-1,1,100); CDFpair = []; CDFnonpair = []; CDFpair_dis = []; CDFnonpair_dis = [];
for ii=1:100,
    CDFpair(ii) = sum([Sp.rho_nd]<=rrr(ii))/Npair;
    CDFnonpair(ii) = sum([Snp.rho_nd]<=rrr(ii))/Nnonpair;
    CDFpair_dis(ii) = sum([Sp_dis.rho_disatt]<=rrr(ii))/Npair_dis;
    CDFnonpair_dis(ii) = sum([Snp_dis.rho_disatt]<=rrr(ii))/Nnonpair_dis;
end
kstest2(CDFnonpair_dis, CDFpair_dis, 0.01), kstest2(CDFnonpair_dis, CDFpair_dis, 0.05),
kstest2(CDFnonpair, CDFpair, 0.01), kstest2(CDFnonpair, CDFpair, 0.05),
plot(rrr, CDFnonpair); xplot(rrr, CDFpair, 'linewidth', 3);
plot(rrr, CDFnonpair_dis); xplot(rrr, CDFpair_dis, 'linewidth', 3);
set(gca,'fontsize', 12)
legend({'non-pairs' 'pairs'}, 'location', 'northwest');
xlabel('Whole-audiogram, nondiagonal xcorr'); ylabel('CDF')
xlim([-0.6 0.6]); ylim([-0.02 1.05]);
text (0.1, 0.7, strvcat('======KS test======', 'CDFs differ at 5% sign. level, ', 'but not at 4% sign. level'));
title('All contralateral responses of cells from pairs');

% === distr of pairs vs nonpairs: detrended, ipsi+contra together
S2nonsame = S2([S2.tt_both_1] & ~[S2.isSame]);
qpair = [S2nonsame.isPair]; qnonpair = ~qpair;
Sp = S2nonsame(qpair); Snp = S2nonsame(qnonpair); Npair = numel(Sp); Nnonpair = numel(Snp);
Sp_dis = Sp([Sp.tt_both_1]); Snp_dis = Snp([Snp.tt_both_1]); Npair_dis = numel(Sp_dis); Nnonpair_dis = numel(Snp_dis);
rrr = linspace(-1.1,1.1,1000); CDFpair = []; CDFnonpair = []; CDFpair_dis = []; CDFnonpair_dis = [];
for ii=1:numel(rrr),
    CDFpair(ii) = sum([Sp.rho_nd]<=rrr(ii))/Npair;
    CDFnonpair(ii) = sum([Snp.rho_nd]<=rrr(ii))/Nnonpair;
    CDFpair_dis(ii) = sum([Sp_dis.rho_disatt]<=rrr(ii))/Npair_dis;
    CDFnonpair_dis(ii) = sum([Snp_dis.rho_disatt]<=rrr(ii))/Nnonpair_dis;
end
kstest2([Sp_dis.rho_disatt], [Snp_dis.rho_disatt], 0.01), kstest2([Sp_dis.rho_disatt], [Snp_dis.rho_disatt], 0.05), 
kstest2([Sp.rho_nd], [Snp.rho_nd], 0.01), kstest2([Sp.rho_nd], [Snp.rho_nd], 0.05), 
plot(rrr, CDFnonpair); xplot(rrr, CDFpair, 'linewidth', 3);  xlim([-0.6 0.6]); 
%plot(rrr, CDFnonpair_dis); xplot(rrr, CDFpair_dis, 'linewidth', 3); xlim([-1 1]); 
set(gca,'fontsize', 12)
legend({'non-pairs' 'pairs'}, 'location', 'northwest');
xlabel('Whole-audiogram, nondiagonal xcorr'); ylabel('CDF')
ylim([-0.02 1.05]);
text (min(xlim)+0.1, 0.6, strvcat('======KS test======', 'CDFs do not differ at 5% sign. level'));
title('All contralateral responses of cells from pairs');


%=======crosscovariance spectrum & crosscor=======
irep = 2;
[Ya, Ncond, dt] = readRueRep('p81213A', irep, 10, 1); [Yb, Ncond, dt] = readRueRep('p81213B', irep, 10, 1);
Ya0 = Ya(:, end-2*29+1:end); Yb0 = Yb(:, end-2*29+1:end); % only the 0-dB SPL condition
% cross power spec 
[Cxx, Freq_c] = cpsd(Ya0(:)-mean(Ya0(:)),Yb0(:)-mean(Yb0(:)),[],[],[], 25e2); xplot(Freq/1e3, P2dB(Cxx), 'k', 'linewidth', 2);
[psdA, Freq] = pwelch(Ya0(:)-mean(Ya0(:)),[],[],[], 25e2); xplot(Freq/1e3, P2dB(psdA), 'r')
[psdB, Freq] = pwelch(Yb0(:)-mean(Yb0(:)),[],[],[], 25e2); xplot(Freq/1e3, P2dB(psdB), 'g')
xlog125; xlabel('Frequency (kHz)'); ylabel('Power (dB)');
% xcorr fcns
set(figure,'units', 'normalized', 'position', [0.12 0.396 0.877 0.504])
[XC, ilag] = xcorr(Ya0(:)-mean(Ya0(:)), Yb0(:)-mean(Yb0(:)), 4e4, 'coeff');
[ACa, ilag] = xcorr(Ya0(:)-mean(Ya0(:)), Ya0(:)-mean(Ya0(:)), 4e4, 'coeff');
[ACb, ilag] = xcorr(Yb0(:)-mean(Yb0(:)), Yb0(:)-mean(Yb0(:)), 4e4, 'coeff');
plot(ilag*dt, XC); xplot(ilag*dt, ACa, 'r'); xplot(ilag*dt, ACb, 'g'); 
legend({'Xcorr' 'autocorr A' 'autocorr B'});
TracePlotInterface(gcf);
% does not look like artefacts (non-neural stuff) play any role in the
% high diagonal correlation of p81213A-p81213B. Most likely is the true
% neural response to ambient acoustic noise, ir some "internally generated"
% common excitation.

% SUPFIG final panel 
RL=sort(RueList); RL = RL(1:2:end);
for ipair=5, 
    S = RuePopcorr(RL{2*ipair-1},RL{2*ipair}); 
    qqn = S.rho; ir = 1:1+size(qqn,1):numel(qqn); qqd = qqn(ir);  qqn(ir)=[]; 
    ah1=subplot(2,1,1); hist(qqn); xlim([-1 1]);
    title([S.FN1 '/' S.FN2]);
    ah2=subplot(2,1,2); hist(qqd); xlim([-1 1]);
%     pause;
%     clf;
end

qq = load('D:\Data\RueData\IBW\Including_Zero_dB\SupFig1_A-F')


%======LARGER DATASET DEC 2010===================================
%=======population stuff=============================
RL = Ruelist; % all recordings
RLc = RL(~cellfun(@isempty, regexp(Ruelist, 'c'))); % only the contra recs
RL2 = cellfun(@(str)str(1:length(str)-1), RLc, 'UniformOutput', false) % remove trailing 'c'
S2 = RueCorrPop(RL2,1); % comparing across all pairs, identicals included. Exclude 0-dB conditions!
save D:\Data\RueData\IBW2\whole_audiogram_xcorr.mat S2
load D:\Data\RueData\IBW2\whole_audiogram_xcorr.mat 
% === KS test on distr of pairs vs nonpairs: detrended, ipsi+contra together
S2nonsame = S2([S2.tt_both_1] & ~[S2.isSame]);
S2nonsame = S2(~[S2.isSame]);
qpair = [S2nonsame.isPair]; qnonpair = ~qpair;
Sp = S2nonsame(qpair); Snp = S2nonsame(qnonpair); Npair = numel(Sp); Nnonpair = numel(Snp);
Sp_dis = Sp([Sp.tt_both_1]); Snp_dis = Snp([Snp.tt_both_1]); Npair_dis = numel(Sp_dis); Nnonpair_dis = numel(Snp_dis);
rrr = linspace(-1.1,1.1,1000); CDFpair = []; CDFnonpair = []; CDFpair_dis = []; CDFnonpair_dis = [];
for ii=1:numel(rrr),
    CDFpair(ii) = sum([Sp.rho_nd]<=rrr(ii))/Npair;
    CDFnonpair(ii) = sum([Snp.rho_nd]<=rrr(ii))/Nnonpair;
    CDFpair_dis(ii) = sum([Sp_dis.rho_disatt]<=rrr(ii))/Npair_dis;
    CDFnonpair_dis(ii) = sum([Snp_dis.rho_disatt]<=rrr(ii))/Nnonpair_dis;
end
KS2_1_p_np = kstest2([Sp.rho_nd], [Snp.rho_nd], 0.01), kstest2([Sp.rho_nd], [Snp.rho_nd], 0.05), 
KS2_1_disatt_p_np = kstest2([Sp_dis.rho_disatt], [Snp_dis.rho_disatt], 0.01), 
KS2_5_disatt_p_np = kstest2([Sp_dis.rho_disatt], [Snp_dis.rho_disatt], 0.05), 
f1; plot(rrr, CDFnonpair); xplot(rrr, CDFpair, 'linewidth', 3);  xlim([-0.6 0.6]); 
set(gca,'fontsize', 12)
legend({'non-pairs' 'pairs'}, 'location', 'northwest');
xlabel('Whole-audiogram, nondiagonal xcorr'); ylabel('CDF')
ylim([-0.02 1.05]);
text (min(xlim)+0.1, 0.6, strvcat('======KS test======', 'CDFs do not differ at 5% sign. level'));
title('All contralateral responses of cells from pairs');
f2; plot(rrr, CDFnonpair_dis); xplot(rrr, CDFpair_dis, 'linewidth', 3); xlim([-1 1]); 
set(gca,'fontsize', 12)
legend({'non-pairs' 'pairs'}, 'location', 'northwest');
xlabel('Whole-audiogram, nondiagonal xcorr'); ylabel('CDF')
ylim([-0.02 1.05]);
text (min(xlim)+0.1, 0.6, strvcat('======KS test======', 'CDFs do not differ at 5% sign. level'));
title('All contralateral responses of cells from pairs');

RL = Ruelist; % all recordings
RLc = RL(~cellfun(@isempty, regexp(Ruelist, 'c'))); % only the contra recs
RLi = RL(~cellfun(@isempty, regexp(Ruelist, 'i'))); % only the ipsi recs
for ii=1:round(numel(RLc)/2),
    RueTriplot(RLc{2*ii-1}, RLc{2*ii});
    pause
    aa;
end
for ii=1:round(numel(RLi)/2),
    RueTriplot(RLi{2*ii-1}, RLi{2*ii});
    pause;
    aa;
end
clear SZ
for ii=1:round(numel(RLi)/2),
    ii
    SZ(ii) = Rue_Zero_dB_Xcorr(RLi{2*ii-1}(1:1:end-1), RLi{2*ii}(1:1:end-1)); % omit trailing 'i'
end
RueXlswrite(SZ, 'D:\Data\RueData\IBW2\zero_dB\zero_dB_Xcorr.xls');
for ii=1:round(numel(RLi)/2),
    RueHexplot(RLi{2*ii-1}(1:1:end-1), RLi{2*ii}(1:1:end-1),1); % omit trailing 'i'
    pause(3); aa;
end

%========microhist
RL = Ruelist; % all recordings
RLc = RL(~cellfun(@isempty, regexp(Ruelist, 'c'))); % only the contra recs
for ii=1:round(numel(RLc)/2),
    FN = RLc{2*ii-1}(1:end-2)
    RueHist(FN);
    pause(1);
    aa;
end
% Rue's selction of microhist
Q = RueHist('p101028'); aa;
clear QQ
BinC = linspace(-1,1,20);
X = Q.XAA(:,:,17); % 9.1 kHz
idiag = 1:(1+size(X,1)):numel(X);
Pnd = X; Pnd(idiag)=[]; Pd = X(idiag);
A_nd = hist(Pnd, BinC);
A_d = hist(Pd, BinC);
QQ.mean_A_nd = mean(Pnd);
QQ.mean_A_d = mean(Pd);
QQ.var_A_nd = var(Pnd);
QQ.var_A_d = var(Pd);
%=
X = Q.XBB(:,:,17); % 9.1 kHz
idiag = 1:(1+size(X,1)):numel(X);
Pnd = X; Pnd(idiag)=[]; Pd = X(idiag);
B_nd = hist(Pnd, BinC);
B_d = hist(Pd, BinC);
QQ.mean_B_nd = mean(Pnd);
QQ.mean_B_d = mean(Pd);
QQ.var_B_nd = var(Pnd);
QQ.var_B_d = var(Pd);
%=
X = Q.XAB(:,:,17); % 9.1 kHz
idiag = 1:(1+size(X,1)):numel(X);
Pnd = X; Pnd(idiag)=[]; Pd = X(idiag);
AB_nd = hist(Pnd, BinC);
AB_d = hist(Pd, BinC);
QQ.mean_AB_nd = mean(Pnd);
QQ.mean_AB_d = mean(Pd);
QQ.var_AB_nd = var(Pnd);
QQ.var_AB_d = var(Pd);
clear TTT
[TTT(1:20).BinCenter] = DealElements(BinC);
[TTT(1:20).A_nd] = DealElements(A_nd);
[TTT(1:20).A_d] = DealElements(A_d);
[TTT(1:20).B_nd] = DealElements(B_nd);
[TTT(1:20).B_d] = DealElements(B_d);
[TTT(1:20).AB_nd] = DealElements(AB_nd);
[TTT(1:20).AB_d] = DealElements(AB_d);
RueXlswrite(TTT, 'D:\Data\RueData\IBW2\microhist\histo_p101028_9kHz_80dB')
RueXlswrite(QQ, 'D:\Data\RueData\IBW2\microhist\mean_var_p101028_9kHz_80dB')







