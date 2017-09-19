%============JL varstat grapes============
S = JLvarStats;
% correct bookkeeping error! see below
Sp = JLfindpartners(S);
Sp = JLgetdate(Sp);
qBp = ([Sp.StimType]=='B' & ~cellfun('isempty', {Sp.I_partners}) & ~cellfun('isempty', {Sp.C_partners}));
[Sp.hasMonPartners] = dealelements(qBp);

save D:\processed_data\JL\JLvarStats\varstats_with_partners Sp


% binInt
iBp = find([Sp.hasMonPartners]);
B = [];
for ii=1:numel(iBp),
    sb = Sp(iBp(ii));
    ui = sb.I_partners;
    uc = sb.C_partners;
    M = min(numel(ui), numel(uc));
    for jj=1:M,
       b = JLbinint(Sp([Sp.UniqueRecordingIndex]==ui(jj)), Sp([Sp.UniqueRecordingIndex]==uc(jj)), sb);
       B = [B, b];
    end
    aa;
end
save D:\processed_data\JL\JLvarStats\bin_int B

% the following loads Sp & B as computed above
load D:\processed_data\JL\JLvarStats\varstats_with_partners.mat
load D:\processed_data\JL\JLvarStats\bin_int B
% general varStat bookkeeping
CellIdx = unique([Sp.UniqueCellIndex]);
iI = find([Sp.StimType]=='I');
iC = find([Sp.StimType]=='C');
iB = find([Sp.StimType]=='B');
iBp = find([Sp.hasMonPartners]);

%---------------------------------------------------
% ----variance accounted for by stimulation --------
%---------------------------------------------------

% ----scatter plot ipsi vs contra  -------
set(figure,'units', 'normalized', 'position', [0.28 0.277 0.641 0.622])
set(gca, 'fontsize', 12);
plot([Sp(iI).ipsi_varAcc], [Sp(iI).contra_varAcc], '.'); 
xplot([Sp(iC).ipsi_varAcc], [Sp(iC).contra_varAcc], 'r.'); 
xplot([Sp(iB).ipsi_varAcc], [Sp(iB).contra_varAcc], 'g.'); 
xplot([Sp(iBp).ipsi_varAcc], [Sp(iBp).contra_varAcc], 'm.'); 
xlim([0 100]); ylim(xlim); grid on; axis square;
xlabel('% var accounted for by Ipsi stim', 'fontsize', 14);
ylabel('% var accounted for by Contra stim', 'fontsize', 14);

%----varAcnt summed between contra and ipsi------
set(figure,'units', 'normalized', 'position', [0.347 0.278 0.608 0.62]);
set(gca, 'fontsize', 12);
Iacc = [Sp(iBp).ipsi_varAcc];
Cacc = [Sp(iBp).contra_varAcc];
Tacc = Cacc+Iacc;
hist(Tacc, 0:5:100);
xlim([0 100]);
Istr = sprintf('Ipsi: %0.1f +/- %d %%', mean(Iacc), round(std(Iacc)));
Cstr = sprintf('Contra: %0.1f +/- %d %%', mean(Cacc), round(std(Cacc)));
Tstr = sprintf('C+I: %0.1f +/- %d %%', mean(Tacc), round(std(Tacc)));
text(0.15, 0.8, Istr, 'units', 'normalized');
text(0.15, 0.75, Cstr, 'units', 'normalized');
text(0.15, 0.7, Tstr, 'units', 'normalized');
xlabel('Contra+Ipsi varAcc (%)', 'fontsize', 14);
ylim([0 170]);

%----varAcnt subtracted between contra and ipsi------
set(figure,'units', 'normalized', 'position', [0.347 0.278 0.608 0.62]);
set(gca, 'fontsize', 12);
Iacc = [Sp(iBp).ipsi_varAcc];
Cacc = [Sp(iBp).contra_varAcc];
Bacc = Cacc-Iacc;
hist(Bacc, -100:5:100);
xlim([-100 100]);
Istr = sprintf('Ipsi: %0.1f +/- %d %%', mean(Iacc), round(std(Iacc)));
Cstr = sprintf('Contra: %0.1f +/- %d %%', mean(Cacc), round(std(Cacc)));
Bstr = sprintf('C-I: %0.1f +/- %d %%', mean(Bacc), round(std(Bacc)));
text(0.15, 0.8, Istr, 'units', 'normalized');
text(0.15, 0.75, Cstr, 'units', 'normalized');
text(0.15, 0.7, Bstr, 'units', 'normalized');
xlabel('Contra-Ipsi varAcc (%)', 'fontsize', 14);
ylim([0 300]);

%----VarAcnt Bin versus I+C
set(figure,'units', 'normalized', 'position', [0.28 0.277 0.641 0.622])
set(gca, 'fontsize', 12);
plot([Sp(iB).ipsi_varAcc] + [Sp(iB).contra_varAcc], [Sp(iB).bin_varAcc], 'k.'); 
xplot([Sp(iBp).ipsi_varAcc] + [Sp(iBp).contra_varAcc], [Sp(iBp).bin_varAcc], 'm.'); 
xlim([0 100]); ylim(xlim); grid on; axis square;
xlabel('% I+C var accounted for', 'fontsize', 14);
ylabel('% Bin var accounted for', 'fontsize', 14);
xplot(xlim,xlim,'b');

%----most varAcnt per cell
set(figure,'units', 'normalized', 'position', [0.28 0.277 0.641 0.622])
set(gca, 'fontsize', 12);
PeakSp = [];
for cellidx = CellIdx,
    sp = Sp([Sp.UniqueCellIndex]==cellidx);
    [dum, imax] = max([sp.bin_varAcc]);
    PeakSp = [PeakSp sp(imax)];
    disp(sp(imax).UniqueRecordingIndex);
end
h = plot([PeakSp.Freq1], [PeakSp.bin_varAcc], 's');
qdisp = @(dum,sp)[sp.ExpID '-' num2str(sp.icell) '  ' num2str(sp.Freq1) ' Hz, ' num2str(sp.SPL1) ' dB  (' num2str(sp.UniqueRecordingIndex) ')'];
IDpoints(h, '', PeakSp, qdisp, 'JLbeatVar', @(dum,sp)JLbeatVar(sp));

%----variation over time of baseline magnitude
iplot = 0;
for cellidx = CellIdx,
    iplot=iplot+1;
    figure(1+floor((iplot-1)/12));
    set(gcf, 'PaperOrientation', 'landscape');
    ipan = 1+rem(iplot-1,12);
    subplot(3,4,ipan);
    sp = Sp([Sp.UniqueCellIndex]==cellidx);
    sp = sortAccord(sp,[sp.abfdatenum]);
    dn = ([sp.abfdatenum] - sp(1).abfdatenum)*(24*60); % minutes since first rec
    plot(dn, [sp.sil_totVar], '*-');
    ylim([0 max(ylim)]);
    if iplot==7, ylim([0 0.1]);
    elseif iplot==10, ylim([0 0.02]);
    elseif iplot==27, ylim([0 0.25]);
    elseif iplot==33, ylim([0 0.03]);
    elseif iplot==39, ylim([0 0.1]);
    elseif iplot==41, ylim([0 0.04]);
    elseif iplot==43, ylim([0 0.02]);
    end
    title([sp(1).ExpID ' cell ' num2str(sp(1).icell)]);
    if ipan==9, xlabel('Time since first recording (minutes)'); end;
    if ipan==5, ylabel('Variance of spont act (mV^2)'); end
end
subsaf

% ----frequency dependency of VarAcc
plot([Sp(iBp).Freq1], [Sp(iBp).ipsi_varAcc]+[Sp(iBp).contra_varAcc], '.');
xlabel('Frequency (Hz)');
ylabel('% VarAcc I+C');
ylim([-0.5 100]);
xlim([0 1600]);

%  ----binInt pop stats
set(figure,'units', 'normalized', 'position', [0.197 0.354 0.718 0.439])
subplot(2,3,1);
hist([B.Rmax_ipsi],100)
xlim([0 1]);
title('ipsi');
subplot(2,3,4);
hist([B.Rmax_contra],100)
xlim([0 1]);
xlabel('Max Corr');
title('contra');
%
taubin = linspace(-1000,1000,80);
subplot(2,3,2);
xlim([-1000 1000]);
hist([B.tau_ipsi],taubin)
xlim([-1000 1000]);
subplot(2,3,5);
hist([B.tau_contra],taubin)
xlim([-1000 1000]);
xlabel('Lag (\mus)');
%
Yi = [B.Yi]; Yc = [B.Yc]; Yb = [B.Yb];
rati = P2dB([Yi.sil_totVar]./[Yb.sil_totVar]);
ratc = P2dB([Yc.sil_totVar]./[Yb.sil_totVar]);
clear Yi Yc Yb
dbbin = -20:2:20;
subplot(2,3,3);
hist(rati, dbbin);
xlim([-20 20]);
subplot(2,3,6);
hist(ratc, dbbin)
xlim([-20 20]);
xlabel('Spont magn mon/bin (dB)');

%---binint: rho vs mag ratio
fisherz = @(rho)log((1+rho)./(1-rho));
subplot(1,2,1);
plot(rati, fisherz([B.Rmax_ipsi]), '.');
xlabel('Spont magn mon/bin (dB)');
ylabel('FisherZ(rho)');
subplot(1,2,2);
plot(ratc, fisherz([B.Rmax_contra]), '.');
xlabel('Spont magn mon/bin (dB)');
ylabel('FisherZ(rho)');

%====================fix--- 204/11 bookkeeping error
JLreadBeats('RG10204',11);
JJJ = [JbB_40dB  JbB_50dB  JbB_60dB  JbC_40dB  JbC_50dB  JbC_60dB  JbI_40dB  JbI_50dB  JbI_60dB];
clear Jb*

ihit = intersect(strmatch('RG10204', {Sp.ExpID}), find([Sp.icell]==11));
AllSeriesIndex = unique([Sp(ihit).UniqueSeriesIndex]);
for ii=1:numel(AllSeriesIndex),
   sidx = AllSeriesIndex(ii);
   iser = find([Sp.UniqueSeriesIndex]==sidx); % indices of this series in Sp
   iJJJ = find([JJJ.UniqueSeriesIndex]==sidx);
   seriesID = JJJ(iJJJ(1)).seriesID, % the correct seriesID
   [Sp(iser).seriesID] = deal(seriesID);
end


% =====JLanova aside=======
MostVarAcc = [178010605  178035401  191010604  191020902  191031402  195034302  ...
    196020611  196041803  197050601  197081904  198010101  198020612  201010102  ...
    201020301  201031507  204106005  204125805  204136905  204147806  209020306  ...
    209030705  209041105  209051705  209061803  214010104  214020502  214030604  ...
    214041302  214051401  214061902  216010404  216021704  219010204  219021202  ...
    219034202  219044504  219054902  219065303  223010302  223020704  223031105  ...
    223041204  223051804  223063001  223073105];
for ii=1:numel(MostVarAcc),
    Spb = Sp([Sp.UniqueRecordingIndex]==MostVarAcc(ii));
    JLanova(Spb);
    pause;
    aa;
    drawnow;
end
% =========new style of addressing data sets=============
D = JLdbase
structview(D([D.iexp]==214 & [D.icond]==1))
% nice example of multi-peakedness, binaural interaction & spike timing
www=JLdatastruct({'RG10214' '6-8' 2})
JLanova(www.I_partners(1))
JLanova(www.C_partners(1))
JLanova(www.B_partners(1))
www=JLdatastruct({'RG10214' '6-8' 2})
qq=JLgetBeatStruct(www)
subplot(2,1,1);
plot(qq.Tsnip, qq.Snip(:,1:10))
subplot(2,1,2);
dSnip = smoothen(diff(qq.Snip), 0.25, qq.dt);
plot(qq.Tsnip(2:end), dSnip(:,1:10)/qq.dt)


%=========align spikes==========
JLspikeAlign({'RG10214' '6-8' 2}, 0.1, [-0.6 -0.125]);

%======attempt at event-based analysis==========
JLevents({'RG10214' '6-8' 2},0.2,1.5,[550 3000], 10);

%=============cell props==============
structview(D([D.iexp]==191 & [D.icond]==1));
Dc = D([D.iseries]==1 & [D.icond]==1); % representative dataset; 1 per cell
for ii=1:45, JDS(ii) = JLdatastruct(Dc(ii)); end
JDS(4).MHcomment = '==';
MHcomm = {JDS.MHcomment};
[strvcat(num2str([Dc.iexp]')) repmat('-',45,1) strvcat(num2str([Dc.icell]')) repmat(' ',45,1) strvcat(MHcomm{:})]
% ans =
% 178- 1 Weakly binaural. Spikes only at 100 Hz.                            
% 178- 3 decent binaural cell                                               
% 191- 1 huge, atypical spikes. Only weakly binaural.                       
% 191- 2 ==                                                                 
% 191- 3 good AP thr                                                        
% 195- 3 well defined thr                                                   
% 196- 2 decent thr                                                         
% 196- 4 mediocre thr definition                                            
% 197- 5 pretty reliable thr definition                                     
% 197- 8 immense spikes; trivial to trigger them.                           
% 198- 1 good AP thr; bookkeeping error: JbB_70dB misses one ABF            
% 198- 2 good AP thr                                                        
% 201- 1 good AP thr                                                        
% 201- 2 good AP thr                                                        
% 201- 3 good AP thr                                                        
% 204-10 AP thr is amazingly stable over recordings                         
% 204-11 poorly defined thr                                                 
% 204-12 reasonably well defined thr                                        
% 204-13 reasonably well defined thr                                        
% 209- 2 reasonably well defined thr                                        
% 209- 3 reasonably well-defined thr; JbB_60dB has huge spikes              
% 209- 4 reasonably well-defined thr                                        
% 209- 5 reasonably well-defined thr                                        
% 209- 6 reasonably well-defined thr                                        
% 214- 1 reasonably well-defined thr                                        
% 214- 2 hard to define a thr w this inverted (?) waveform                  
% 214- 3 not a sharp thr                                                    
% 214- 4 not a sharp thr                                                    
% 214- 5 strange waveforms: sinks                                           
% 214- 6 okayish                                                            
% 216- 1 reasonable thr def                                                 
% 216- 2 decent thr def                                                     
% 219- 1 decent thr def                                                     
% 219- 2 poorly def thr; hardly spiking, but binaural waveforms             
% 219- 3 weak binaurality; mostly ipsi driven; very few spikes. well def thr
% 219- 4 reasonable thr                                                     
% 219- 5 reasonable thr; good S/N                                           
% 219- 6 okay thr                                                           
% 223- 1 need high thr to correct quantization noise                        
% 223- 2 okay thr                                                           
% 223- 3 okay thr                                                           
% 223- 4 not a very crisp thr                                               
% 223- 5 not a very crisp thr                                               
% 223- 6 not a very crisp thr                                               
% 223- 7 not a well defined thr                                             













%=====interpretation of countour plots=========
ph=linspace(0,2*pi);
X = cos(ph); Y = cos(ph).'; 
contourf(Y*X)
aa
contourf(ph,ph,Y*X)
aa
contourf(ph/2/p,ph/2/pi,Y*X)
contourf(ph/2/pi,ph/2/pi,Y*X)
aaa
[X,Y] = SameSize(X,Y);
contourf(ph/2/pi,ph/2/pi, X.*Y)
aa
contourf(ph/2/pi,ph/2/pi, X+Y)
aa
contourf(ph/2/pi,ph/2/pi, (X+Y).^abs(X+Y))
contourf(ph/2/pi,ph/2/pi, (X+Y)./abs(X+Y))
contourf(ph/2/pi,ph/2/pi, (X+Y)./abs(X+Y).^0.2)
Z = (X+Y)./abs(X+Y).^0.2;
Zx = mean(Z,2); Zy = mean(Z,1); [Zx,Zy] = SameSize(Zx,Zy);
contourf(ph/2/pi,ph/2/pi, Z)
contourf(ph/2/pi,ph/2/pi, Zx+Zy)
aa
contourf(ph/2/pi,ph/2/pi, Z-(Zx+Zy))




