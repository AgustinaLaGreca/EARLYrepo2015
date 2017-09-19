function [Fr, MG, PH, Lnoise, W, isortFreq] = ZWOAEsubspec(D, idataset, subType, Ref);
% ZWOAEsubspec - partial spectrum of ZWOAE data
%   [Freq, Magn, Phase, Lnoise,W, isortF] = ZWOAEspec(iGerbil, idataset, SubType, Ref) 
%   returns the spectrum of ZWOAE data (iGerbil, idataset) of a particular 
%   subgroup of the recorded ZWOAE spectrum. 
%   SubType is one of (quotes omitted):
%       zwuis: zwuis components in stimulus.
%      zwuis2: 2nd harmonics of zwuis components.
%      zwuis3: 3rd harmonics of zwuis components.
%      single: the single primary in the stimulus
%  suppressor: the suppressor tone [if present].
%    mimicker: the mimicker tone [if present].
%        near: components of near group (see ZWOAEmatrices)
%         far: components of far group  (see ZWOAEmatrices)
%       suplo: components of lower suppression group (see ZWOAEmatrices)
%       suphi: components of higher suppression group (see ZWOAEmatrices)
%      supall: all components of suppression group (see ZWOAEmatrices)       
%         all: all DP components (see ZWOAEmatrices)
%  Ref is the reference used for evaluating the phases. Possible values
%  for Ref are:
%     'stim': [default] Reference phases are taken from the stimulus 
%             specification ("electrical stimulus phases").
%     'resp': Reference phases are taken from the data. Note that this will
%             result in all phases being zero for the SubType='zwuis' case.
%          X: if an numerical array X is specified, X(k) will be used as a
%             reference for the component at frequency Freq(k). The phases
%             in X must be specifed in cycles.
%  In the two former cases, the phases from the stimulus components are
%  used to compute the reference phases of the DPs ("2*phi1-ph2" etc). In
%  the last case, the user-specified reference phases correspond directly
%  to the DP components.
%  
%  Output args:
%        Fr: sorted (!) frequency [kHz] of requested components
%        MG: magnitude (dB) of requested components
%        PH: phase (cycles) of requested components
%    Lnoise: estimated noise level (dB) at requested components
%         W: weight factors for fitting spectral components. W is computed
%            from the MG and Lnoise using the function SNR2W.
%    isortF: index array "unsorting" all previous arrays to make them
%            compatible with ZWOAEmatrices (if applicable). That is,
%            Fr(isortF) are the frequencies in the same order as M*Fstim, 
%            where M is the appropriate matrix returned by ZWOAEmatrices.
%
%  If there is only one output arg, S, then all the above outputs are
%  collected in a single struct, and many more .. just try it! In this case
%  (single output arg S), you may also specify multiple-component idataset,
%  in which case a struct array is returned.
% 
%  See also getZWOAEdata, addZWOAEsubspec, ZWOAEfit, SNR2W.

if nargin<4, Ref = 'stim'; end

if ischar(Ref),
    [Ref, Mess] = keywordMatch(Ref, {'stimulus' 'response'});
elseif ~isnumeric(Ref), error('Ref input arg must be valid keyword or numerical array.');
end

if numel(idataset)>1 && nargout<2, % use recursion to get struct array
    for ii=1:numel(idataset),
        S(ii) = ZWOAEsubspec(D, idataset(ii), subType, Ref);
    end
    Fr = reshape(S,size(idataset));
    return;
end


D = ZWOAEimport(D, idataset);
% compute matrix M that computes requested freqs from stimulus freq
Fstim = [D.Fzwuis(:); D.Fsingle; D.Fsup; D.Fmim]; % all stimulus freqs 
[subType, Mess] = keywordMatch(subType, {'zwuis' 'zwuis2' 'zwuis3' 'single' ...
    'suppressor' 'near' 'far' 'suplo' 'suphi' 'supall' 'all' 'mimicker'});
error(Mess);
switch subType,
    case 'single',
        M = [zeros(1,D.Nzwuis) 1 0 0];
    case 'suppressor',
        M = [zeros(1,D.Nzwuis+1) 1 0];
    case 'mimicker',
        M = [zeros(1,D.Nzwuis+1) 0 1];
    case 'zwuis',
        M = [diag(ones(1,D.Nzwuis)) zeros(D.Nzwuis,3)];
    case 'zwuis2',
        M = 2*[diag(ones(1,D.Nzwuis)) zeros(D.Nzwuis,3)];
    case 'zwuis3',
        M = 3*[diag(ones(1,D.Nzwuis)) zeros(D.Nzwuis,3)];
    otherwise,
        M = getFieldOrDefault(ZWOAEmatrices(D.Nzwuis),['M' subType], []);
        M = [M zeros(size(M,1),2)]; % zero cols to accomodate fsup & fmim
end
% M,Fstim
[Fr, isortF] = sort(M*Fstim);
[MG, icomp] = SpecComp(D.df, D.MG, Fr);
PH = SpecComp(D.df, D.PH, Fr);
% take reference phases into account
if isnumeric(Ref), % user specified - accept as is
    PHref = Ref(:);
else, % take from stimulus or response
    switch Ref, % get reference phases of the stimulus components (see help text)
        case 'stimulus',
            PHrefStim = [D.PHzwuis(:); D.PHsingle; D.PHsup; D.PHmim];
            if isempty(PHrefStim),
                error('Cannot use stimulus phases as a reference; they are not stored in dataset.');
            end
        case 'response',
            PHrefStim = SpecComp(D.df, D.PH, Fstim);
    end
    PHref = M*PHrefStim; % stim phase -> DP phase
    PHref = PHref(isortF); % same sorting as applied earlier to Fr
end
PH = PH-PHref;
PM = phaseMean(PH);
PH = (PH-PM)+0.5; PH = mod(PH,1)-0.5+PM;
Pnoise = D.Pnf;
%Lnoise = polyval(Pnoise,Fr);
Lnoise = ZWOAEgetNoise(D, Fr);
W = SNR2W(MG-Lnoise,0,20); % restrict dynamic range of weighting factors
W = SNR2W(MG-Lnoise);
isortFreq(isortF) = 1:max(isortF); % isortFreq inverse is permutation of isortF

if nargout<2, % collect outputs in struct
    MGsummed = P2dB(sum(dB2P(MG))); % total power of requested components
    ExpID = D.ExpID; recID = D.recID; RecType = D.RecType; 
    Fzwuis = D.Fzwuis; Fz_mean = D.Fz_mean; Fsingle = D.Fsingle; Fsup = D.Fsup; Fmim = D.Fmim;
    Lzwuis = D.Lzwuis; Lsingle = D.Lsingle; Lsup = D.Lsup; Lmim = D.Lmim;
    Fr = CollectInStruct(ExpID, recID, RecType, '-', ...
        Fzwuis, Fz_mean, Fsingle, Fsup, Fmim, Lzwuis, Lsingle, Lsup, Lmim, '-', ...
        subType, Ref, Fr, MG, PH, Lnoise, W, isortFreq, MGsummed, icomp);
end
    
