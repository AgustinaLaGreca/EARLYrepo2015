function PR = JLpredictMonRate(Uidx_B, Uidx_I, Uidx_C);
% JLpredictMonRate - predict monaural spike rates from binaural recs
%    P = JLpredictMonRate(Uidx_B, Uidx_I, Uidx_C)
%  OR
%    P = JLpredictMonRate(CSC), with
%    CSC a struct containing fields Uidx_B, Uidx_I, Uidx_C, e.g.,
%    an element of
%    retrieve_dbase('D:\processed_data\JL\JLsplitvar\JLsplitvar.dbase')
%
% The output of JLpredictMonRate is stored in database
%  D:\processed_data\JL\PredMonSpikes\JLpredictMonRate_XXX.dbase,
%  where XXX stands for the BinCount parameter (see code).
%  Retrieve it by, e.g.,
%  DBFN = 'D:\processed_data\JL\PredMonSpikes\JLpredictMonRate_10.dbase';
%  PR = retrieve_dbase(DBFN);


BinCount = 10; % min # spikes in each V bin

if isstruct(Uidx_B),
    S = Uidx_B;
    [Uidx_B, Uidx_I, Uidx_C] = deal(S.Uidx_B, S.Uidx_I, S.Uidx_C);
    clear S;
end
UniqueTripletIndex = -Uidx_B;

DBFN = ['D:\processed_data\JL\PredMonSpikes\JLpredictMonRate_' num2str(BinCount) '.dbase'];
init_dbase(DBFN, 'UniqueTripletIndex', 'onlyifnew');


PLB = JLpowerlaw(Uidx_B);
PLI = JLpowerlaw(Uidx_I);
PLC = JLpowerlaw(Uidx_C);
% true number of spikes
[Nsp_B, Nsp_I, Nsp_C] = deal(PLB.NspikeStim, PLI.NspikeStim, PLC.NspikeStim);

[Vbin, Vedge, pSpike] = local_repower(PLB, BinCount);

% input var & peak-peak
[Var_B, Var_I, Var_C] = local_CSfield('VarDriv', Uidx_B, Uidx_I, Uidx_C);
[SpontVar_B, SpontVar_I, SpontVar_C] = local_CSfield('VarSpont1', Uidx_B, Uidx_I, Uidx_C);
[P2P_Imon, P2P_Ibin] = local_CSfield('P2PMeanIpsiCycle', Uidx_I, Uidx_B);
[P2P_Cmon, P2P_Cbin] = local_CSfield('P2PMeanContraCycle', Uidx_C, Uidx_B);

% vector strengths
[VS_Imon, VS_Ibin] = local_SSfield('VS_I', Uidx_I, Uidx_B);
[VS_Cmon, VS_Cbin] = local_SSfield('VS_C', Uidx_C, Uidx_B);
VS_B = local_SSfield('VS_B', Uidx_B);
[VSalpha_Imon, VSalpha_Ibin] = local_SSfield('vs_alpha_I', Uidx_I, Uidx_B);
[VSalpha_Cmon, VSalpha_Cbin] = local_SSfield('vs_alpha_C', Uidx_C, Uidx_B);
VSalpha_B = local_SSfield('vs_alpha_B', Uidx_B);


% plot(Vbin, pSpike, '*-')
% xplot(PLB.Vinst, PLB.cond_spike_prob, 'rs-');

Nb = local_predictNsp(PLB, PLB.MeanDriv, Vbin-PLB.MeanDriv, pSpike);
scaleFactor = Nsp_B/Nb; % normalize "autoprediction"
Nsp_pred_I = round(scaleFactor*local_predictNsp(PLI, PLI.MeanDriv, Vbin-PLB.MeanDriv, pSpike));
Nsp_pred_C = round(scaleFactor*local_predictNsp(PLC, PLC.MeanDriv, Vbin-PLB.MeanDriv, pSpike));

id = local_ID(Uidx_B);

pr1 = CollectInStruct(UniqueTripletIndex, Uidx_B, Uidx_I, Uidx_C);
pr2 = CollectInStruct('-Params', BinCount, ...
    '-SpikeCounts', Nsp_B, Nsp_I, Nsp_C, ...
    '-Predictions', Nsp_pred_I, Nsp_pred_C, scaleFactor, ...
    '-InputVar', Var_B, Var_I, Var_C, ...
    '-SpontVar', SpontVar_B, SpontVar_I, SpontVar_C, ...
    '-P2P', P2P_Imon, P2P_Cmon, P2P_Ibin, P2P_Cbin, ...
    '-VS', VS_Imon, VS_Cmon, VS_Ibin, VS_Cbin, VS_B, ...
    '-VSalpha', VSalpha_Imon, VSalpha_Cmon, VSalpha_Ibin, VSalpha_Cbin, VSalpha_B ...
    );

PR = structJoin(pr1, id, pr2);
Add2dbase(DBFN, PR);

%==============================================
function [Vbin, Vedge, pSpike] = local_repower(P, BinCount);
% redo powerlaw computation, using a balanced bin occupancy
if P.NspikeStim<BinCount, % senseless exercise
    [Vbin, Vedge, pSpike] = deal(nan);
    return;
end
VspB = sort(P.mean_at_spt); % V at spike times
VspB = sort(P.mean_at_spt); % V at spike times
Nsp = P.NspikeStim;
group_edge = BinCount*(0:floor(Nsp/BinCount));
group_edge = group_edge - max(group_edge) + Nsp;
group_edge(1) = 0; % now the borders demarkate equal amount of spikes, except the first one, which can have more
group_edge = 0.5+group_edge;
Vedge = interp1(1:Nsp, VspB, group_edge);
[VminSp VmaxSp] = minmax(P.mean_at_spt);
meanDist = mean(diff(denan(Vedge)));
Vedge(end) = VmaxSp+0.01*meanDist;
Vedge(1) = VminSp-0.01*meanDist;
% prepend bins that cover the nonspiking-domain
Vmin = min(P.mean2d(:))-0.01*meanDist;
Nstep = ceil((Vedge(1)-Vmin)/meanDist);
if Nstep>1,
    Vextra = linspace(Vmin, Vedge(1), Nstep);
    Vedge = [Vextra, Vedge(2:end)];
end
Vedge(end) = max(P.mean2d(:))+0.01*meanDist;
hist_mean_at_spt = histc(VspB, Vedge)
hist_mean_at_spt(end) = [];
histV = histc(P.mean2d(:).', Vedge);
histV(end) = [];
pSpike = hist_mean_at_spt./(P.DurPerCell*histV);
Vbin = (Vedge(1:end-1)+Vedge(2:end))/2;

function N = local_predictNsp(P, Vbaseline, Vbin, pSpike);
% predict # spikes from "power law" pSpike vs Vbin
if isnan(Vbin),
    N = nan;
else,
    Vinst = P.mean2d(:).'-Vbaseline;
    p_inst = interp1(Vbin, pSpike, Vinst, 'linear', 'extrap');
    N = sum(p_inst)*P.DurPerCell;
end

function id = local_ID(Uidx);
DS = JLdatastruct(Uidx);
FNS = {'ID_______________________' 'iexp' 'icell' 'icell_run' 'iseq' 'SPL' 'chan' 'iseries' 'iseries_run' 'icond' 'freq'};
id = structpart(DS, FNS);

function varargout = local_CSfield(FN, varargin);
for ii=1:numel(varargin),
    CS = JLcycleStats(varargin{ii});
    varargout{ii} = CS.(FN);
end

function varargout = local_SSfield(FN, varargin);
for ii=1:numel(varargin),
    CS = JLspikeStats(varargin{ii});
    varargout{ii} = CS.(FN);
end










