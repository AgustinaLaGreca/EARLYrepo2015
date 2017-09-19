function [DP, MG] = supspec(DS, Chan, iCond);
% dataset/subspec - spectral analysis of SUP data
%    supspec(DS, Chan, iCond) produces a magnitude spectrum of the SUP data in DS. 
%    Chan is the channel spec of the analog channel. See dataset/anachan.
%    iCond is the stimulus condition (counting)
%
%    [DP, MG] = supspec(...) returns info on the DPs in struct array DP and
%    on the magnitude spectrum in MG.
%  
%    See also dataset/Magn_Phase_Plot.

Rcrit = 0.001; % Rayleighcriterion for significance of phase locking

if ~isequal('SUP', upper(stimtype(DS))),
    error('Dataset is not a SUP dataset.');
end

[D dt]= anamean(DS,Chan);
NsamCyc = DS.Stim.NsamZwuisPeriod; % Number of samples per cycle

% Correct D for offset samples from TDT
NsamSpec = NsamCyc*floor(size(D,1)/NsamCyc);
D = D(1:NsamSpec,:);

% Compute mean of zwuis repetitions
qq = simplegate(D(:,iCond),round(NsamCyc/sqrt(2))); % apply temporal window, but avoid using stimulus periodicity
qq = qq - mean(qq(:));
[qq M] = LoopMean(qq,NsamCyc);

% Compute the FFT
Spec = fft(qq);
Magn = A2dB(abs(Spec));
Phase = cangle(Spec);
df = M/(dt*size(D,1)); % frequency spacing in kHz
Freq = Xaxis(Spec,df);
Alpha = anarayleigh(qq, dt, 1e3*Freq, 0, 10);

% auxiliary stuff
Fsup = DS.Stim.Fsup(:)/1e3; % suppressor freq in kHz
Fprobe = DS.Stim.Fprobe(iCond)/1e3; % probe freq in kHz
Fprim = [Fsup(:); Fprobe]; % all stim freqs ("primary frequencies") in col array
isam = @(freq)1+round(freq/df); % freq (kHz) -> sample index
Nsup = numel(Fsup);
PhiPrimStim = [DS.Stim.SupStartPhase(:); 0]; % electrical phases of primaries (to be used as references)
PhiPrimResp = Phase(isam(Fprim)); % phases of linear responses to primaries (to be used as references)

DP = []; % will become struct array whose elements describe a DP class
%ClassNames = words2cell('Sum/Diff/SumSup/DiffSup/Supn/Upper/Lower/Sup3/Diff2/Sum2',  '/');
ClassNames = Words2cell('Sum/Diff/SumSup/DiffSup/Supn/Upper/Sup3/Diff2/Sum2',  '/');
for idp=1:numel(ClassNames),
    try,
        dp.class = ClassNames{idp};
        dp.M = local_DPmatrix(Nsup, dp.class);
        dp.Freq = dp.M*Fprim; % frequencies of DP family
        dp.isam = isam(dp.Freq);
        dp.Alpha = Alpha(dp.isam);
        dp.Magn = Magn(dp.isam);
        dp.PhaseRefStim = dp.M*PhiPrimStim; % phase references from response
        dp.PhaseRefResp = dp.M*PhiPrimResp; % phase references from response
        dp.rawPhase = Phase(dp.isam);
        dp.PhaseReStim = ucunwrap(dp.rawPhase-dp.PhaseRefStim, dp.Freq,'-tozero');
        dp.PhaseReResp = ucunwrap(dp.rawPhase-dp.PhaseRefResp, dp.Freq,'-tozero');
        dp.qsig = dp.Alpha<=Rcrit;
        [freq, phi] = sortAccord(dp.Freq(dp.qsig), dp.PhaseReResp(dp.qsig), dp.Freq(dp.qsig));
        if sum(dp.qsig)>1,
            P = polyfit(freq, phi,1);
            dp.tauResp = -P(1)/1e3; % group delay in us
            dp.ph0Resp = P(2);
        else,
            dp.tauResp = nan;
            dp.ph0Resp = nan;
        end
        DP = [DP, dp];
    catch
        warning(dp.class);
    end % try/catch
end


% =========plot========
FplotMax = max(2*Fprobe)+3.5*mean(Fsup);
PS = local_plotStruct;
STL = stimlist(DS);
set(figure,'units', 'normalized', 'position', [0.0313 0.222 0.93 0.678]);
Str = strvcat([expname(DS) ' ' num2str(irec(DS)) ' ' IDstring(DS)], ...
    [STL.Spec2 ';  CF = ' num2str(DS.Stim.SupCenterFreq/1e3) ' kHz; BW = ' num2str(DS.Stim.SupBandWidth) ' Hz' ] ,...
    [STL.SPL ';  CF = ' num2str(round(Fprobe)) ' kHz'] ); 

%--magn--
MG = CollectInStruct(df, Magn, isam, Fsup, Fprobe);
ha(1) = subplot(2,1,1);
dplot(df, Magn, PS.Spec);
xplot(Fsup, Magn(isam(Fsup)), PS.Supr);
xplot(Fprobe, Magn(isam(Fprobe)), PS.Probe);
for idp=1:numel(DP),
    dp = DP(idp);
    xplot(dp.Freq, dp.Magn+pmask(dp.Alpha<=Rcrit), PS.(dp.class));
end
xlim([0 FplotMax]);
ylim(max(Magn(2:end))+[-100 10]);
text(0.5, 0.9, Str, 'units', 'normaliz');
%--phase--
ha(2) = subplot(2,1,2);
for idp=1:numel(DP),
    local_phiplot(DP(idp), 'Resp', PS, Rcrit);
end
set(ha(2), 'ytick', -10:0.25:10);
ylim([-0.75 0.75]);
grid on
linkaxes(ha,'x')


%=======================================================

function local_phiplot(dp, Ref, PS, Rcrit);
fn = ['PhaseRe' Ref]; % DP name
Phi = dp.(fn)+pmask(dp.Alpha<=Rcrit);
xplot(dp.Freq, -1+Phi, PS.(dp.class));
xplot(dp.Freq, Phi, PS.(dp.class));
xplot(dp.Freq, 1+Phi, PS.(dp.class));
[freq, phi] = sortAccord(dp.Freq(dp.qsig), Phi(dp.qsig), dp.Freq(dp.qsig));
if numel(freq>2),
    P = polyfit(freq, phi,1);
    xplot(freq, polyval(P,freq), 'k-');
    text(mean(freq), mean(phi)+0.2, [num2str(round(-1e3*P(1))) ' \mus'], 'horizontalalign', 'center');
end

function PS = local_plotStruct;
% spectrum
PS.Spec = struct('color',[1 1 1]*.7,'LineWidth',2);
% noise floor
PS.Floor = struct('color',[1 1 1]*.2,'LineWidth',2);
% probe
PS.Probe = struct('color','r', 'markerfacecolor','r', 'marker', 'o',...
    'linestyle', 'none', 'markersize', 4, 'linewidth',0.5);
% suppressor
PS.Supr = struct('color','r', 'markerfacecolor','w', 'marker', 'o',...
    'linestyle', 'none', 'markersize', 4, 'linewidth',0.5);
% fprobe + Fsup
PS.Sum = struct('color','m', 'markerfacecolor','m', 'marker', 'p',...
    'linestyle', 'none', 'linewidth',0.5);
% % fprobe - Fsup
PS.Diff = struct('color',[0 0 0.7], 'markerfacecolor',[0 0 0.7],...
    'marker', 'p', 'linestyle', 'none', 'linewidth',0.5);
% fprobe + Fsup
PS.SumSup = struct('color','m', 'markerfacecolor','w', 'marker', 'p',...
    'linestyle', 'none', 'linewidth',0.5);
% % fprobe - Fsup
PS.DiffSup = struct('color',[0 0 0.7], 'markerfacecolor','w',...
    'marker', 'p', 'linestyle', 'none', 'linewidth',0.5);
% suppression: Fprobe + diffmatrix(Fsup)
PS.Supn = struct('color','g', 'markerfacecolor','g', 'marker', 's',...
    'linestyle', 'none', 'markersize', 4, 'linewidth',0.5);
% upper: Fprobe + Fsup(i) + Fsup(j)
PS.Upper = struct('color','b', 'markerfacecolor','b', 'marker', '^',...
    'linestyle', 'none', 'markersize', 4, 'linewidth',0.5);
% lower: Fprobe - Fsup(i) - Fsup(j)
PS.Lower = struct('color','b', 'markerfacecolor','b', 'marker', 'v',...
    'linestyle', 'none', 'markersize', 4, 'linewidth',0.5);
% sup3: Fsup(i)+Fsup(j)-Fsup(k)
PS.Sup3 = struct('color','g', 'markerfacecolor','w', 'marker', 's',...
    'linestyle', 'none', 'markersize', 4, 'linewidth',0.5);
% sum2: 2*Fprobe+Fsup(i)
PS.Sum2 = struct('color','m', 'markerfacecolor','m', 'marker', '>',...
    'linestyle', 'none', 'markersize', 4, 'linewidth',0.5);
% diff2: 2*Fprobe-Fsup(i)
PS.Diff2 = struct('color',[0 0 0.7], 'markerfacecolor',[0 0 0.7], 'marker', '<',...
    'linestyle', 'none', 'markersize', 4, 'linewidth',0.5);


function M = local_DPmatrix(Nsup, nam);
% matrix to be multiplied to col array [Fsup(:); Fprobe] yielding named group
% of DP frequencies. Fsup in ascending order.
switch lower(nam),
    case 'diffsup', % fsup(i) - fsup(j)
        M = local_Xcol(-DPmatrix(Nsup,2,0),0); % this yields pos Fsup differences du to ascending Fsup order
    case 'sumsup',  % fsup(i) + fsup(j)
        M = local_Xcol(DPmatrix(Nsup,2,2),0); 
    case 'sup3',    % fsup(i) + fsup(j)-fsup(k)
        M = local_Xcol(DPmatrix(Nsup,3,1),0); 
    case 'diff',    % fprobe  - fsup(i)
        M = local_Xcol(-eye(Nsup),1); 
    case 'sum',     % fprobe + fsup(i)
        M = local_Xcol(eye(Nsup),2); 
    case 'diff2',    % fprobe  - fsup(i)
        M = local_Xcol(-eye(Nsup),2); 
    case 'sum2',     % fprobe + fsup(i)
        M = local_Xcol(eye(Nsup),1); 
    case 'supn',    % fprobe + fsup(i) - fsup(j)
        M = local_Xcol([-DPmatrix(Nsup,2,0); DPmatrix(Nsup,2,0)],1);
    case 'lower',     % fprobe - fsup(i) - fsup(j)
        M = local_Xcol(-DPmatrix(Nsup,2,2),1);
    case 'upper',     % fprobe + fsup(i) + fsup(j)
        M = local_Xcol(DPmatrix(Nsup,2,2),1);
    otherwise,
        error(['Unkown DP class ''' nam '''.']);
end

function M = local_Xcol(M, val);
% append one column with constant value
M = [M, repmat(val, size(M,1),1)];







