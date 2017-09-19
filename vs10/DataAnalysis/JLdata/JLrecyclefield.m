function JLrecyclefield(Uidx);
% JLrecyclefield - receptive-field view of cycle waveforms
%    JLrecyclefield(Uidx)

% identify datasets
db = JLdbase(Uidx);
DB = JLdbase;
DB = DB([DB.icell_run]==db.icell_run); % select data from same cell 
% inventory of all stim conditions
freq = unique([DB.freq]);
SPL = unique([DB.SPL]);
Nfreq = numel(freq);
Nspl = numel(SPL);
XX = (0:Nfreq-1)/Nfreq; XX = 0.05+0.9*XX;
DX = 0.8*mean(diff(XX));
if isnan(DX), DX = 0.9; end
%YY = (0:Nspl-1)/Nspl; YY = 0.05+0.85*YY;
YY = (0:Nspl-1)/Nspl; YY = 0.08+0.82*YY;
DY = 0.75*mean(diff(YY));
if isnan(DY), DY = 0.9; end

if Nspl==1, FigHeight = 0.2; Yfig = 0.7;
elseif Nspl==2, FigHeight = 0.4; Yfig = 0.55;
elseif Nspl==3, FigHeight = 0.6; Yfig = 0.35;
else, FigHeight = 0.8; Yfig = 0.15;
end
FigPos = [0.0141 Yfig 0.978 FigHeight];
set(figure,'units', 'normalized', 'position', FigPos);
Ymin = inf; Ymax = -inf;
for ispl=Nspl:-1:1,
    for ifreq=1:Nfreq,
        Apos = [XX(ifreq) YY(ispl) DX DY];
        ha(ispl,ifreq) = axes('position', Apos);
        LegStr = {};
        db = DB([DB.freq]==freq(ifreq) & [DB.SPL]==SPL(ispl));
        db_i = db([db.chan]=='I');
        db_c = db([db.chan]=='C');
        db_b = db([db.chan]=='B');
        for ii=1:numel(db_i),
            W = JLwaveforms(db_i(ii).UniqueRecordingIndex);
            xdplot(W.dt_IpsiMean, W.IpsiMeanrec-mean(W.IpsiMeanrec), 'b');
            LegStr = [LegStr 'I/I'];
        end
        for ii=1:numel(db_c),
            W = JLwaveforms(db_c(ii).UniqueRecordingIndex);
            xdplot(W.dt_ContraMean, W.ContraMeanrec-mean(W.ContraMeanrec), 'r');
            LegStr = [LegStr 'C/C'];
        end
        for ii=1:numel(db_b),
            W = JLwaveforms(db_b(ii).UniqueRecordingIndex);
            xdplot(W.dt_IpsiMean, W.IpsiMeanrec-mean(W.IpsiMeanrec), 'b', 'linewidth', 2);
            LegStr = [LegStr 'I/B'];
            xdplot(W.dt_ContraMean, W.ContraMeanrec-mean(W.ContraMeanrec), 'r', 'linewidth', 2);
            LegStr = [LegStr 'C/B'];
        end
        xlim([0, 1e3/freq(ifreq)]);
        Ymin = min(min(ylim), Ymin);
        Ymax = max(max(ylim), Ymax);
        set(gca, 'xticklabel', ''); 
        if ifreq>1, set(gca, 'yticklabel', ''); end
        if ifreq==1,
            text(0.5, 0.9, sprintf('%d dB',SPL(ispl)), 'horizontalalign', 'center', 'units', 'normalized');
        end
        if ispl==1,
            text(0.5, -0.1, sprintf('%d Hz', round(freq(ifreq))), 'horizontalalign', 'center', 'units', 'normalized');
        end
        if ifreq==round(Nfreq/2) && ispl==Nspl,
            title(sprintf('exp %d cell %d', DB(1).iexp, DB(1).icell), 'fontsize', 14);
        end
    end
end
set(ha(:),'ylim', [Ymin Ymax]);









