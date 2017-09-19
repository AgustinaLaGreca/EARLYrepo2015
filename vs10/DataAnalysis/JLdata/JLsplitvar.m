function S = JLsplitvar(Uidx, doPlot, varargin);
% JLsplitvar - splitting the variance of a (I,C,B) triplet of recordings
%   S = JLsplitvar([Uidx_I Uidx_C Uidx_B], doPlot)
%   S = JLsplitvar([Uidx_I Uidx_C Uidx_B], doPlot, ...plotargs...)
%
%   S = JLsplitvar(Uidx_B) finds the monaural partners and does the whole
%   series, i.e. all freqs shared by BI Nad MON recs.
%
%  Example: JLsplitvar([214062402, 214062302, 214062202]);
%
% yyy = @(Uidx) sprintf('%% JLsplitvar(%d,1, ''marker'', ploma(1), lico(1)); %% %s;', Uidx, strrep(JLtitle(Uidx), ' 100 Hz ',''))
%
% JLsplitvar(195034301,1, 'marker', ploma(1), lico(1)); % Exp 195 cell 3 --- B   65 dB  (series 1); ***
% JLsplitvar(195034401,1, 'marker', ploma(2), lico(2)); % Exp 195 cell 3 --- B   75 dB  (series 2); **
%
% JLsplitvar(196020502,1, 'marker', ploma(1), lico(1)); % Exp 196 cell 2  --- B   80 dB  (series 1);
%
% JLsplitvar(198021702,1, 'marker', ploma(1), lico(1)); % Exp 198 cell 2  --- B   40 dB  (series 4);
% JLsplitvar(198020402,1, 'marker', ploma(2), lico(2)); % Exp 198 cell 2  --- B   70 dB  (series 2); **
%
% JLsplitvar(201020702,1, 'marker', ploma(1), lico(1)); % Exp 201 cell 2  --- B  150 Hz   40 dB  (series 2);
% JLsplitvar(201020402,1, 'marker', ploma(2), lico(2)); % Exp 201 cell 2 --- B  150 Hz   60 dB  (series 3); ***
%
% JLsplitvar(204112501,1, 'marker', ploma(1), lico(1)); % Exp 204 cell 10 --- B   50 dB  (series 8);
% JLsplitvar(204112801,1, 'marker', ploma(2), lico(2)); % Exp 204 cell 10  --- B   60 dB  (series 2);
% JLsplitvar(204113101,1, 'marker', ploma(3), lico(3)); % Exp 204 cell 10 --- B   70 dB  (series 10);
%
% JLsplitvar(204126401,1, 'marker', ploma(1), lico(1)); % Exp 204 cell 11  --- B   40 dB  (series 1);
% JLsplitvar(204126301,1, 'marker', ploma(2), lico(2)); % Exp 204 cell 11  --- B   50 dB  (series 2);
% JLsplitvar(204125801,1, 'marker', ploma(3), lico(3)); % Exp 204 cell 11 --- B   60 dB  (series 3); **
%
% JLsplitvar(204137401,1, 'marker', ploma(1), lico(1)); % Exp 204 cell 12  --- B   40 dB  (series 2);
% JLsplitvar(204137101,1, 'marker', ploma(2), lico(2)); % Exp 204 cell 12  --- B   50 dB  (series 3);
% JLsplitvar(204136701,1, 'marker', ploma(3), lico(3)); % Exp 204 cell 12 --- B   60 dB  (series 4); ***
%
% JLsplitvar(204148601,1, 'marker', ploma(1), lico(1)); % Exp 204 cell 13  --- B   50 dB  (series 2);
% JLsplitvar(204147801,1, 'marker', ploma(2), lico(2)); % Exp 204 cell 13  --- B   60 dB  (series 3);
%
% JLsplitvar(209020301,1, 'marker', ploma(1), lico(1)); % Exp 209 cell 2 --- B   60 dB  (series 1); ***
%
% JLsplitvar(209041401,1, 'marker', ploma(1), lico(1)); % Exp 209 cell 4  --- B   60 dB  (series 1)
%
% JLsplitvar(214041001,1, 'marker', ploma(1), lico(1)); % Exp 214 cell 4 --- B   60 dB  (series 2); **
%
% JLsplitvar(214063701,1, 'marker', ploma(1), lico(1)); % Exp 214 cell 6  --- B   10 dB  (series 3)
% JLsplitvar(214063401,1, 'marker', ploma(2), lico(2)); % Exp 214 cell 6  --- B   20 dB  (series 4)
% JLsplitvar(214063101,1, 'marker', ploma(3), lico(3)); % Exp 214 cell 6  --- B   30 dB  (series 5)
% JLsplitvar(214062801,1, 'marker', ploma(4), lico(4)); % Exp 214 cell 6  --- B   40 dB  (series 6)
% JLsplitvar(214062501,1, 'marker', ploma(5), lico(5)); % Exp 214 cell 6  --- B   50 dB  (series 7)
% JLsplitvar(214061801,1, 'marker', ploma(6), lico(6)); % Exp 214 cell 6 --- B   60 dB  (series 1) ***
% JLsplitvar(214061901,1, 'marker', ploma(7), lico(7)); % Exp 214 cell 6  --- B   70 dB  (series 9)
% JLsplitvar(214062201,1, 'marker', ploma(8), lico(8)); % Exp 214 cell 6  --- B   80 dB  (series 10)
%
% JLsplitvar(216010101,1, 'marker', ploma(1), lico(1)); % Exp 216 cell 1  --- B   70 dB  (series 2)
%
% JLsplitvar(216021401,1, 'marker', ploma(1), lico(1)); % Exp 216 cell 2  --- B   40 dB  (series 5)
% JLsplitvar(216021101,1, 'marker', ploma(2), lico(2)); % Exp 216 cell 2  --- B   50 dB  (series 6)
% JLsplitvar(216020801,1, 'marker', ploma(3), lico(3)); % Exp 216 cell 2  --- B   60 dB  (series 7)
% JLsplitvar(216020501,1, 'marker', ploma(4), lico(4)); % Exp 216 cell 2  --- B   70 dB  (series 8)
%
% JLsplitvar(219031902,1, 'marker', ploma(1), lico(1)); % Exp 219 cell 3  --- B  200 Hz   60 dB  (series 4)
% JLsplitvar(219034103,1, 'marker', ploma(2), lico(2)); % Exp 219 cell 3  --- B  300 Hz   80 dB  (series 6)
%
% JLsplitvar(223041202,1, 'marker', ploma(1), lico(1)); % Exp 223 cell 4  --- B  200 Hz   50 dB  (series% 3)
%

[doPlot] = arginDefaults('doPlot',0);

if numel(Uidx)==1, % find monaural partners; do whole series
    D = JLdatastruct(Uidx);
    DB = JLdbase('iseries_run', D.iseries_run);
    NB = numel(DB);
    S = [];
    for ii=1:NB,
        try,
            D = JLdatastruct(DB(ii).UniqueRecordingIndex);
            Uidx_I = D.I_partners(1);
            Uidx_C = D.C_partners(1);
            if ~isempty(Uidx_I) && ~isempty(Uidx_C),
                s = JLsplitvar([Uidx_I, Uidx_C, DB(ii).UniqueRecordingIndex]);
                S = [S, s];
            end
        catch,
            warning(lasterr);
        end
    end
    if doPlot,
        local_Plot(S, doPlot, varargin{:});
    end
    return
end

[Uidx_I, Uidx_C, Uidx_B] = DealElements(Uidx);

D = JLdatastruct(Uidx_B);
CFN = ['D:\processed_data\JL\JLsplitvar\JLsplitvar_' num2str(D.UniqueSeriesIndex)];
DBFN = 'D:\processed_data\JL\JLsplitvar\JLsplitvar.dbase';
init_dbase(DBFN, 'UniqueTripletIndex', 'onlyifnew');
[S, CFN, CP] = getcache(CFN, [Uidx_I, Uidx_C, Uidx_B]);


if isempty(S), % compute S
   UniqueTripletIndex = -Uidx_B; % by convention, triplets are tagged by minus the BIN Uidx
   iexp = D.iexp;
   icell = D.icell;
   icell_run = D.icell_run;
   iseries = D.iseries;
   iseries_run = D.iseries_run;
   freq = D.freq;
   SPL = D.SPL;
   % ipsi
   W = JLwaveforms(Uidx_I); 
   CSi = JLcycleStats(Uidx_I);
   [minIvar, iminIvar] = min(W.IpsiVar);
   VbaseI = W.IpsiMeanrec(iminIvar);
   RMSmonIpsi = std(W.IpsiMeanrec);
   % contra
   W = JLwaveforms(Uidx_C); 
   CSc = JLcycleStats(Uidx_C);
   [minCvar, iminCvar] = min(W.ContraVar);
   VbaseC = W.ContraMeanrec(iminCvar);
   RMSmonContra = std(W.ContraMeanrec);
   % bin
   W = JLwaveforms(Uidx_B); 
   CSb = JLcycleStats(Uidx_B);
   minBinVar = min(W.BeatVar);
   [SpontVarIpsi, SpontVarContra, SpontVarBin] = deal(CSi.VarSpont1, CSc.VarSpont1, CSb.VarSpont1);
   MeanSpontVar = mean([CSi.VarSpont1 CSc.VarSpont1]);
   RMSbinIpsi = std(W.IpsiMeanrec);
   RMSbinContra = std(W.ContraMeanrec);
   % collect
   S = CollectInStruct(UniqueTripletIndex, Uidx_I, Uidx_C, Uidx_B, ...
       '-ID', iexp, icell, icell_run, iseries, iseries_run, freq, SPL, ...
       '-JLsplitvar', ...
       MeanSpontVar, SpontVarIpsi, SpontVarContra, SpontVarBin, minBinVar, minIvar, iminIvar, ...
       VbaseI, minCvar, iminCvar, VbaseC,...
       '-RMSmeanCycle', RMSmonIpsi, RMSmonContra, RMSbinIpsi, RMSbinContra);
   % cache & dbase
   Add2dbase(DBFN, S);
   putcache(CFN, 50, CP, S);
end

function local_Plot(S, doPlot, varargin);
if ~doPlot, return; end
h = gcf;
if isempty(getGUIdata(gcf, 'FigureType', '')),
    set(h,'tag', mfilename);
    set(h,'units', 'normalized', 'position', [0.28 0.353 0.523 0.547]);
    setGUIdata(h,'FigureType',mfilename);
    setGUIdata(h,'LegendStr',{});
    setGUIdata(h,'LegendAxh',[]);
end
ah = xplot([S.freq], [S.minIvar]./[S.SpontVarIpsi] + [S.minCvar]./[S.SpontVarContra] - [S.minBinVar]./[S.SpontVarBin], 'linestyle', '-', varargin{:});
xplot([S.freq], [S.minBinVar]./[S.SpontVarBin], '.-', varargin{:}, 'markersize', 5);
xplot(xlim, [1 1],'k');
xplot(xlim, [2 2],'k:');
title(JLtitle(S(1).Uidx_B));
LegStr = getGUIdata(h,'LegendStr');
Axh = getGUIdata(h,'LegendAxh');
LegStr = [LegStr [num2str(S(1).SPL) ' dB']];
Axh = [Axh ah];
legend(Axh, LegStr, 'location', 'east');
setGUIdata(h,'LegendStr',LegStr);
setGUIdata(h,'LegendAxh',Axh);
xlabel('Frequency (Hz)');
ylabel('Normalized Variance');



