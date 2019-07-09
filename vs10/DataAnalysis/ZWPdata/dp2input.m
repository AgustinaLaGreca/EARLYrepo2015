function [D, Z] = dp2input(d, Cdelay)
% dp2input - IHC input from DP2 spectrum  of spike train evoked by zwuis stim
%  [D, Z] = dp2input(d, Cdelay)
%        d: dataset (ZW or ZWP stimulus).
%   Cdelay: ms delay compensation (needed for phase unwrapping)

AlphaCrit = 0.001; % max p value for Rayleigh test
NminOccur = 3; % minimum occurrence of primary as parent of DP2s

doPlot = nargout<1;

Cdelay = arginDefaults('Cdelay',3); 

Z = zwuisana(d.Stim, digichan(d,1)); 
St = d.Stim; 
SPT = spiketimes(d);
Nsp = sum(Z.R);% # spikes
f2i = @(f)1+round(1e3*f/Z.df); % kHz => index in spectrum


VScrit = rayleigh_rmin(Nsp, AlphaCrit);

Fprim = St.Fzwuis.'/1e3; % kHz primary freq
iprim = f2i(Fprim);
Nprim = numel(Fprim);

DP2 = DPfreqs(Fprim, 2, 0) ; % difference tones only
Fdiff = [DP2.freq].'; % diff freq
idiff = f2i(Fdiff);

VSprim = Z.Zsp(iprim); % complex vector strength (including phase)
MagPrim = a2db(abs(VSprim));
PhasePrim = cunwrap(cangle(VSprim) - St.StartPhase.' + Cdelay*Fprim);
PMprim = pmask(abs(VSprim)>=VScrit(1));

VSdiff = Z.Zsp(idiff); % complex vector strength (including phase)
MagDiff = a2db(abs(VSdiff));
% remove nonsign components
qsig = abs(VSdiff)>=VScrit(end);
Ndp2 = sum(qsig);
DP2 = DP2(qsig);
Fdiff = [DP2.freq].'; % diff freq
idiff = f2i(Fdiff);
[VSdiff, MagDiff] = deal(VSdiff(qsig), MagDiff(qsig));
Mw = cat(1,DP2.weight); % weight matrix
% phases
PhaseRef = Mw*St.StartPhase.';
PhaseDiff = cangle(VSdiff) - PhaseRef + Cdelay*Fdiff;


% solve MagInp from MagDiff = abs(Mw)*MagInp
% MagInp1 = abs(Mw)\MagDiff;
% MagInp1 = MagInp1 - nanmean(MagInp1);
[MagInp, FitMag] = local_leftdiv(abs(Mw), MagDiff, NminOccur, false, Fdiff, [1 2 3 4 5]);
MagInp = MagInp - nanmean(MagInp);
1;

% solve PhaseInp from PhaseDiff = Mw*PhaseInp
% PhaseInp = [Mw; ones(1, Nprim)]\[PhaseDiff; 0]; % addition "data" forcing sum(PhaseInp)=0, lifting rank deficiency
PhaseInp = [];



ID = d.ID;
D = CollectInStruct(ID, '-param', AlphaCrit, VScrit, Nsp, Cdelay, ...
    '-apple', Nprim, Fprim, VSprim, MagPrim, PhasePrim, PMprim, ...
    Fdiff, VSdiff, MagDiff, PhaseDiff, PhaseRef, ...
    '-dp2apple', Ndp2, Mw, MagInp, FitMag, PhaseInp);

if doPlot,
    local_plot(D);
end


%============================================
function [X, F] = local_leftdiv(Mw, Data, NminOccur, forceZeroMean, Fpoly, PolyPowers);
% solves X and C from Data = M*X + polyval(C, Fpoly)
Ndata = numel(Data);
Nprim = size(Mw,2);
q_not_rare = sum(abs(Mw)) >= NminOccur; % primaries that occur a suffient # times in Mw
M = Mw(:,q_not_rare);
Zfreq = zscore(Fpoly);
if isempty(PolyPowers), % make VDM  construction below will work
    PolyPowers = zeros(1,0);
end
if forceZeroMean,  % add constraint fixing sum(X)
    M = [M; ones(1, size(M,2))]; 
    Data = [Data; 0];
end
VDM = bsxfun(@power, Zfreq, PolyPowers); % VanderMonde matrix for polynomial fit
Npoly = size(VDM,2);
XX = [M, VDM]\Data;
NN = size(M,2); % number of primary comps retrieved
[XX, C] = deal(XX(1:NN), XX(NN+1:end));
X = nan(1,Nprim,1);;
X(q_not_rare) = XX;
F = CollectInStruct(VDM, C, Zfreq, NminOccur, forceZeroMean, Fpoly, PolyPowers);


function local_plot(D);
fh = figure;
set(fh,'units', 'normalized', 'position', [0.0974 0.312 0.28 0.565]);
subplot(2,1,1); 
plot(D.Fdiff, D.MagDiff, 'o'); 
FilterTrend = D.FitMag.VDM*D.FitMag.C;
FilterTrend = FilterTrend-nanmean(FilterTrend)+ nanmean(D.MagDiff);
[dpfreq, FilterTrend] = sortAccord(D.Fdiff, FilterTrend, D.Fdiff);
xplot(dpfreq, FilterTrend, '-', 'markersize',4);

xlabel('Envelope frequency (kHz)');
ylabel('Magnitude (dB re VS=1)');
subplot(2,1,2); 
plot(D.Fprim, D.MagInp, '*')
xlabel('Primary Frequency (kHz)');
ylabel('Magnitude (dB)');


















