function [S, imRec]=poolZWOAEdata(iGerbil, idataset, ignoreCompanions);
% poolZWOAEdata - combine multiple ZWOAErecordings into one virtual dataset
%    S = poolZWOAEdata(iGerbil, iRec) combined the recorded spectra
%    specified in array iRec into a virtual dataset having the same format
%    as the struct returned by a standard call to ZWOAEimport.
%    This amounts to an weighted averaging of the indivual scomplex
%    spectra, with the priviso that all stimulus-related components are
%    phase-aligned before being summed. Stimulus-reated components are
%    any components that can be extracted by ZWOAEsubspec, e.g.,
%    'far', 'near', 'zwuis', 'mimicker', etc.
%
%    [S, ImRec] = poolZWOAEdata(iGerbil, iRec) also returns an imaginary 
%    index ImRec that allows the pooled data to be addressed as if they
%    were really collected in a single recording. They may be retrieved
%    by subsequently calling 
%
%            ZWOAEimport(iGerbil, Imrec);
%
%    Such saved pooled data are available during the rest of the current 
%    MatLab session.
%
%    See also ZWOAEimport, addZWOAsubspec

if nargin<3, ignoreCompanions=0; end
subSpecTypes = {'zwuis' 'zwuis2' 'zwuis3' 'single' 'suppressor' 'mimicker' 'all'};
stimComps = {'zwuis' 'single' 'sup' 'mim'}; % somponents whose phase (PHzwuis, etc) aer set to zero

Cspec = 0; % initialize complex spectrum
TotWeight = 0; icount=0;
for irec=idataset(:).',
    icount = icount+1;
    D = getZWOAEdata(iGerbil, irec);
    Weight = D.Stimparam.recdur; 
    [csp, RecType, companionID(icount)] = local_corrected_Cspec(iGerbil,irec,subSpecTypes);
    isAcoust(icount) = isequal('ACOUSTIC',RecType);
    Cspec = Cspec + Weight*csp;
    TotWeight = TotWeight + Weight;
end
Cspec = Cspec/TotWeight;
if ~ignoreCompanions && (all(isAcoust) || all(~isAcoust)),
    [dum, imRecComp] = poolZWOAEdata(iGerbil,companionID,1); % 1: ignore companions of companions
else, % no clear companion group
    imRecComp = nan;
end

% convert complex spec back to magn & phase, and place combined spectrum in
% first S specified
S = ZWOAEimport(iGerbil,idataset(1));
[S.MG, S.PH] = Cspec2Magn_Phase(Cspec);
S.recID = idataset; % for the record ;)
S.companionID = imRecComp; 
% S.Pnf = -100;

freq = 1e3*(0:S.Nsam-1).'*S.df; % freq of spectral components in Hz
MG = S.MG;
if strcmpi(S.RecType,'acoustic') && isfield(S,'probename') && ~strcmpi(S.probename,'none'),
    ProbeGain = ProbeLoss(S.probename, freq);
    MG = MG + ProbeGain; % plus sign: we "add" the probe transfer
end
S.Pnf = fitNoise(MG, S.Nsam/S.periodicity, S.df);

% set stimulus phases to zero
for ii=1:numel(stimComps),
    fn = ['PH' stimComps{ii}]; % 'PHzwuis' etc
    S.(fn) = 0*S.(fn); % same size, but zero-valued
end

if nargout>1, % save as imaginary dataset
    imRec = imagZWOAErecording(S);
    % tell the iminary companion pooled data their companion ...
    if ~isnan(imRecComp), imagZWOAErecording(0, imRecComp, imRec); end
end

% ======================
function [Cspec, RecType, companionID]=local_corrected_Cspec(iGerbil,irec,subSpecTypes);
% complex spectrum with special phases corrected based on stimulus phases
S = ZWOAEimport(iGerbil, irec);
for ii=1:numel(subSpecTypes), % extract phase-correced components using ZWOAEsubspec & place them back in overall spec
    SubS = ZWOAEsubspec(iGerbil, irec, subSpecTypes{ii}, 'stim');
    S.MG = setSpecComp(S.df, S.MG, SubS.Fr, SubS.MG);
    S.PH = setSpecComp(S.df, S.PH, SubS.Fr, SubS.PH);
end
Cspec = Magn_phase2Cspec(S);
RecType = S.RecType;
companionID = S.companionID;







