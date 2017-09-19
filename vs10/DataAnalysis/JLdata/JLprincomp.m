function [P, isel] = JLprincomp(S, doPlot, iselect, PlotCrit);
% JLprincomp - principal component analysis of JLbeat data
%   JLprincomp(S), where S is output of JLbeat or JLbeatStats.


if nargin<2, doPlot=(nargout<1); end % only plot if no return arg is requested
if nargin<3, iselect=[]; end
if nargin<4, PlotCrit=[]; end

if isempty(iselect), iselect={[] []}; end
if isempty(PlotCrit), PlotCrit=[0 1]; end


[JB, S] = JLgetBeatStruct(S); % make sure to have both "raw" output JB of JLbeat & JLbeatStats output ST
stimOnset = 500; % ms start of stim re start of rec

[D, DS, L]=readTKABF(S);
dt = 1/S.Fsam; % sample period in ms
rec = D.AD(1).samples; % recording
D.AD = []; % save memory
Time = timeaxis(rec,dt); % tim eaxis in ms
measID = [S.ExpID ' ' S.seriesID '(' num2str(S.icond) '):  ' num2str(S.Freq1) ' Hz, ['  trimspace(num2str([S.SPL1 S.SPL2])) '] dB'];
% left
C1 = local_cycleSort(dt,rec, S.Period1);
C2 = local_cycleSort(dt,rec, S.Period2);

C1 = local_reject(C1, stimOnset+S.Ttrans, JB.SPTraw, S.APwindow_start, S.APwindow_end, iselect{1});
C2 = local_reject(C2, stimOnset+S.Ttrans, JB.SPTraw, S.APwindow_start, S.APwindow_end, iselect{2});

[C1.pCF, C1.pSC] = princomp(C1.rec);
[C2.pCF, C2.pSC] = princomp(C2.rec);

P = CollectInStruct(measID, '-', C1, C2);
isel = {[] []};
if doPlot, % display
    % ---basic PCA results
    Nplot = 7;
    fh1 = figure;
    set(fh1,'units', 'normalized', 'position', [0.0484 0.435 0.443 0.51]);
    subplot(4,3,[1 6]);
    dplot(P.C1.dt, P.C1.pSC(:,1:Nplot)); 
    title(JB.TTT);
    xlabel('Time (ms)');
    legend(Words2cell(num2str(1:Nplot)), 'location', 'south', 'fontsize', 6);
    fh2 = figure;
    set(fh2,'units', 'normalized', 'position', [0.56 0.434 0.443 0.51]);
    subplot(4,3,[1 6]);
    dplot(P.C2.dt, P.C2.pSC(:,1:Nplot)); 
    xlabel('Time (ms)');
    % ---histograms of PC weights
    figure(fh1);
    P = local_histplot(P, 1);
    figure(fh2);
    P = local_histplot(P, 2);
    % ---overlay of reconstructed waveforms 
    M=15; 
    set(figure,'units', 'normalized', 'position', [0.05 0.0684 0.438 0.325]); 
    [P, iselect_I] = local_overlay(P, 1, 200, PlotCrit);
    set(figure,'units', 'normalized', 'position', [0.56 0.0664 0.438 0.325]); 
    [P, iselect_C] = local_overlay(P, 2, 200, PlotCrit);
    disp(['% ' P.measID]);
    isel = {iselect_I iselect_C};
end
P.plotcluster = @(P,iear, j1,j2)plot(P.(['C' num2str(iear)]).pCF(:,j1), P.(['C' num2str(iear)]).pCF(:,j2), '.');
% Period1
% 
% AnaDur
% S.Ttrans


%===============================
function C = local_cycleSort(dt, rec, StimCycle);
% chop recording in cycles, each containing integer # samples
%
% first change sample rate to make Period contain integer # samples
Nsam = numel(rec);
NsamCycle = round(StimCycle/dt); % rounded # samples per stimulus cycle
dt_new = dt*(StimCycle/dt)/NsamCycle; % sample period resulting in NsamCycle samples/stimCycle
Time_old = dt*(0:Nsam-1).';
Nsam_new = 1+floor(max(Time_old)/dt_new);
Nsam_new = NsamCycle*floor(Nsam_new/NsamCycle);  %integer # cycles
Time_new = dt_new*(0:Nsam_new-1).';
rec = interp1(Time_old, rec(:), Time_new);
dt = dt_new;
Nsam = Nsam_new;
Ncycle = Nsam/NsamCycle;
rec = reshape(rec,NsamCycle, Ncycle);
t0 = (0:Ncycle-1)*StimCycle; % starting times of cycles
t1 = (1:Ncycle)*StimCycle; % end times of cycles
C = CollectInStruct(dt, rec, Nsam, StimCycle, t0, t1, NsamCycle, Ncycle);

function C = local_reject(C, InitialInterval, SPT, APwindow_start, APwindow_end, iselect);
iok = ([C.t0]>InitialInterval);  % cycles occurring after initial segment
iok = iok(:).';
if ~isempty(SPT),
    APstart = SPT + APwindow_start;
    APend = SPT + APwindow_end;
    % find out which cycles contain (part) of an AP
    histEdges = [C.t0, C.t1(end)];
    Nstart = histc(APstart, histEdges); Nstart = Nstart(:).';
    Nend = histc(APend, histEdges); Nend = Nend(:).';
    %dsize(histEdges, iok, Nstart, Nend)
    iok = iok & (Nstart(1:end-1)==0) & (Nend(1:end-1)==0);
end
% apply selection to fields of C
C.rec = C.rec(:,iok);
C.t0 = C.t0(:,iok);
C.t1 = C.t1(:,iok);
% posthoc selection according to iselect
if ~isempty(iselect),
    C.rec = C.rec(:,iselect);
    C.t0 = C.t0(:,iselect);
    C.t1 = C.t1(:,iselect);
end


function P = local_histplot(P, iear);
Cfld = ['C' num2str(iear)]; % 'C1' or 'C2'
c = P.(Cfld);
for ii=1:20, 
    doPlot = (ii<=6);
    if doPlot,
        subplot(4,3,6+ii);
        hist(c.pCF(:,ii),50);
        hp = findobj(gca,'type', 'patch');
    end
    P.(Cfld).lillie(ii) = lillietest(c.pCF(:,ii), 0.001);
    P.(Cfld).ttest(ii) = ttest(c.pCF(:,ii),0,0.001);
    if doPlot,
        if P.(Cfld).ttest(ii),
            set(hp, 'edgecolor', 'b')
        else,
            set(hp, 'edgecolor', 0.5*[1 1 1])
        end
        if P.(Cfld).lillie(ii),
            set(hp, 'facecolor', 'b');
        else,
            set(hp, 'facecolor', 'w');
        end
    end
end;

function [P, iselect] = local_overlay(P, iear, M, plotcrit);
if nargin<3, 
    M = 0;
end
if nargin<4, 
    plotcrit=1; % plot all
end
Cfld = ['C' num2str(iear)]; % 'C1' or 'C2'
c = P.(Cfld);
if isequal(0,M),
    M = find(~c.lillie & ~c.ttest, 1, 'first');
    if isempty(M), M = numel(c.ttest); end
end
predy=c.pSC(:,1:M)*c.pCF(:,1:M).';
vp=var(predy);
if numel(plotcrit)==1,
    if plotcrit>0,%plotcrit=0.3 means Discard lower 30%
        plotcrit = [plotcrit 1];  
    else, %plotcrit=-0.3 means Select lower 30%
        plotcrit = [0 -plotcrit];  
    end
end
[dum, isort] = sort(vp, 'ascend');
i0 = max(1,round(plotcrit(1)*numel(vp)));
i1 = round(plotcrit(2)*numel(vp));
iselect = isort(i0:i1); % indices of events meeting plotcrit
dplot(c.dt, predy(:, iselect));
title([P.measID ' crit=' strrep(trimspace(num2str(plotcrit)), ' ', '-')], 'Interpreter', 'none')





