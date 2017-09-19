function BM = JLbinintMetrix(I,C,B);
% JLbinintMetrix - metrics derived from JLbinint
%   M=JLbinint(iI, iC, iB), where iX are recording indices of I,C,B.
%   M is struct array w metrics as fields. 
%
%   See also JLbeatBinint.

JB = JLbinint(I,C,B,0); % 0: don't plot

m_ii = JB.Yi.ipsi_MeanWave;
m_cc = JB.Yc.contra_MeanWave;
m_bi = JB.Yb.ipsi_MeanWave;
m_bc = JB.Yb.contra_MeanWave;
m_bb = JB.Yb.bin_MeanWave;
%
v_ii = JB.Yi.ipsi_VarWave;
v_cc = JB.Yc.contra_VarWave;
v_bi = JB.Yb.ipsi_VarWave;
v_bc = JB.Yb.contra_VarWave;
v_bb = JB.Yb.bin_VarWave;

BM.iexp = JB.Yi.iexp;
BM.icell = JB.Yi.icell;
BM.icell_run = JB.Yi.icell_run;
BM.SPL = JB.Yi.SPL;
BM.freq = JB.Yi.Freq1;
BM.iseqi = JB.Yi.iseq;
BM.iseqc = JB.Yc.iseq;
BM.iseqb = JB.Yb.iseq;
BM.Ui = I;
BM.Uc = C;
BM.Ub = B;
BM.plot = @()JLbinint(I,C,B,1);
BM.binint = @()JLbinint(I,C,B,0);
%
BM = local_sep(BM, 'mean_waveform');
BM.ipsi_P2P_ipsi = localP2P(m_ii);   % ipsi mean wv, peak to peak
BM.contra_P2P_contra = localP2P(m_cc); % contra mean wv, peak to peak
BM.bin_P2P_ipsi = localP2P(m_bi);   % ipsi mean wv, peak to peak
BM.bin_P2P_contra = localP2P(m_bc); % contra mean wv, peak to peak
BM.bin_P2P_bin = localP2P(m_bb); % contra mean wv, peak to peak
BM.ipsi_mean_ipsi = mean(m_ii);   
BM.contra_mean_contra = mean(m_cc); 
BM.bin_mean_ipsi = mean(m_bi);  
BM.bin_mean_contra = mean(m_bc); 
BM.bin_mean_bin = mean(m_bb); 
%
BM = local_sep(BM, 'variance');
BM.ipsi_totvar = JB.Yi.ipsi_totVar;
BM.ipsi_var_not_ipsi = JB.Yi.ipsi_resVar;
BM.ipsi_spont_var = JB.Yi.sil_totVar;
BM.contra_totvar = JB.Yc.contra_totVar;
BM.contra_var_not_contra = JB.Yc.contra_resVar;
BM.contra_spont_var = JB.Yc.sil_totVar;
BM.bin_totvar = JB.Yb.bin_totVar;
BM.bin_var_not_ipsi = JB.Yb.ipsi_resVar;
BM.bin_var_not_contra = JB.Yb.contra_resVar;
BM.bin_var_not_stim = JB.Yb.bin_resVar;
BM.bin_spont_var = JB.Yb.sil_totVar;
%
BM = local_sep(BM, 'mon_vs_bin');
BM.Rmax_ipsi = JB.Rmax_ipsi;
BM.tau_ipsi = JB.tau_ipsi;
BM.Rmax_contra = JB.Rmax_contra;
BM.tau_contra = JB.tau_ipsi;
%
BM = local_sep(BM, 'peaks');
BM = local_peaks(BM, 'ipsi', 'ipsi', JB.Yi.ipsi_dt, m_ii);
BM = local_peaks(BM, 'contra', 'contra', JB.Yc.contra_dt, m_cc);
BM = local_peaks(BM, 'bin', 'ipsi', JB.Yb.ipsi_dt, m_bi);
BM = local_peaks(BM, 'bin', 'contra', JB.Yb.contra_dt, m_bc);

%======================================
function p2p = localP2P(Y)
p2p = max(Y)-min(Y);

function B = local_sep(B, Str);
B.(['x______' Str '____']) = '____________';

function BM = local_peaks(BM, recEar, AnaEar, dt, Y);
NpFld = [recEar '_Npeak_' AnaEar];
VpFld = [recEar '_Vpeak_' AnaEar];
TpFld = [recEar '_Tpeak_' AnaEar];
TmxFld = [recEar '_Tmax_' AnaEar];
[Tp, Yp] = localmax(dt, Y, 0.4);
BM.(NpFld) = numel(Tp);
BM.(TpFld) = Tp;
BM.(VpFld) = Yp - mean(Y);
[dum, imax] = max(Y);
BM.(TmxFld) = dt*(imax-1);












