function [D, DD, LrefD, LrefDD] = ZWOAEdominance(MG,df,Fzw,Fs,DPtype);
% ZWOAEdominance - spectral dominance metric of DP group
%   [D, DD] = ZWOAEdominance(MG,df,Fzw,Fs,DPtype) computes the dominance 
%    of a given DP type in their spectral region.
%     Inputs
%        MG,df: magn spectrum in dB and freq spacing in kHz (see ZWOAEspec)
%       Fzw,Fs: stimulus freq in kHz (Fzw is zwuis group; Fs single tone)
%       DPtype: type of DP considered. One of 'near' 'far' 'suplo' 'suphi'
%               'supall' 'all'
%    Outputs
%      D: S/N ratio in dB, where "signal power" is total power of the DPtype
%         components and "noise power" is total power the M strongest
%         components in the same spectral region with the exception of any 
%         stimulus components. M is the number of DPtype components.
%     DD: same as D, but now other known DPs ('all') are also excluded from
%         the noise power.
% 
%   [D, DD, LrefD, LrefDD] = ZWOAEdominance(MG,df,Fzw,Fs,DPtype) also
%   returns the average level [dB] of the M strongest "noise" components
%   and the DD version of it (see above).
%
%    ZWOAEdominance(D, DPtype)
%       where D is the struct returned by ZWOAEimport, is the same as
%       ZWOAEdominance(D.MG, D.df, D.Fzwuis, D.Fsingle, DPtype)
%
%
%    See also getZWOAEdata, ZWOAEimport, ZWOAEspec.

if nargin==2, % extract params from data struct (see help text)
    [MG,df,Fzw,Fs,DPtype] = deal(MG.MG, MG.df, MG.Fzwuis, MG.Fsingle, df);
end

freq = df*(0:numel(MG)-1); % freq axis of spectrum
[freq, Fzw, Fs] = local_roundFreq(df, freq, Fzw, Fs);
Fprim = [Fzw(:); Fs]; % primary freqs in col array
Nzwuis = numel(Fzw);

MM = ZWOAEmatrices(Nzwuis);
% select DPtype spectral components & "all" DPs
MMdp = MM.(['M' DPtype]);
MMdpAll = MM.Mall;
Fdp = local_roundFreq(df,MMdp*Fprim);
Fdp = abs(Fdp);
FdpAll = local_roundFreq(df,MMdpAll*Fprim);
MGdp = SpecComp(df, MG, Fdp);
% select reference ("noise") spectral components
Ndp = numel(MGdp); % # DP cmpts
[FdpMin, FdpMax] = minmax(Fdp);
Fref = freq((freq>FdpMin)&(freq<FdpMax)); % all freqs in range ..
FrefD = setdiff(Fref, [Fprim; Fdp]); % .. except primary freqs 
FrefDD = setdiff(Fref, [Fprim; FdpAll]); % .. except primary & DP freq fre
MGrefD = SpecComp(df,MG,FrefD); 
MGrefDD = SpecComp(df,MG,FrefDD); 
% the powers that be
Pdp = sum(dB2P(MGdp)); % the total power of DP compts
PrefD = local_maxpower(MGrefD,Ndp); % power of Ndp strongest ref components
PrefDD = local_maxpower(MGrefDD,Ndp); % ditto, DD version (see help text)
D = 0.1*round(10*P2dB(Pdp/PrefD));
DD = 0.1*round(10*P2dB(Pdp/PrefDD));
LrefD = P2dB(PrefD/Ndp); % mean power per strong ref component
LrefDD = P2dB(PrefDD/Ndp); % ditto, DD version


% ================
function varargout=local_roundFreq(df, varargin);
% round freqs for disambiguation
for ii=1:nargin-1,
    varargout{ii} = 0.01*df*round(100*varargin{ii}/df);
end

function P=local_maxpower(MG,N);
% total power of N strongest components or equivalent
M = numel(MG);
P = sort(dB2P(MG)); % power of components in ascending order
if M==0,
    P=nan;
elseif M<N,
    P = P2dB(N/M)*sum(P); % take total power, correct for shortage of cmpnts
else,
    P = sum(P(end-N+1:end)); % summed power of N strongest components
end




