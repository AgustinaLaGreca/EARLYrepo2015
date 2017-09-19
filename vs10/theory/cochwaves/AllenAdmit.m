function [Y, S] = AllenAdmit(freq, S_BM, fr_TM, S_TM, eta, fr_STS, S_HB);
%  AllenAdm - admittance of section of Allen's 1980 cochlea.
%    Y = AllenAdmit(freq, fr_TM, fr_STS, alpha, beta)
%    Inputs
%        freq: frequency (Hz)
%        S_BM: BM stiffness per longitudinal length unit (Pa)
%       fr_TM: TM resonance freq (Hz)
%        S_TM: TM stiffness per longitudinal length unit (Pa)
%         eta: damping coeff of TM
%      fr_STS: cutoff frequency (Hz) of subtectorial space & OHC stiffness
%        S_HB: OHC hair bundle stiffness per unit length (Pa)

iomega = 2*pi*i*freq;
[om_TM, om_STS] = deal(2*pi*fr_TM, 2*pi*fr_STS);  % angular freqs

Zsts = S_HB./iomega.*(1+iomega./om_STS);
Ztm = S_TM./iomega.*(1+(iomega./om_TM).^2 + eta.*iomega./om_TM);

Zoc = impara(Zsts, Ztm);
Z = Zoc + S_BM./iomega;
Y = 1./Z;
S = CollectInStruct(Zsts, Ztm, Zoc);
