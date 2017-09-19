function SH = JLspontHist(Uidx, doPlot)
%  JLspontHist - plot history of spont activity of givn MSO cell
%     SH = JLspontHist(Uidx, doPlot)


doPlot = arginDefaults('doPlot', nargout<1);

DB = JLdbase;
db = DB([DB.UniqueRecordingIndex]==Uidx);
DB = DB([DB.icell_run]==db.icell_run); % all recs of same cell

[SH, CFN, CP] = getcache(fullfile(processed_datadir, 'JL', 'JLpop', 'sponthist'), db.icell_run);
if isempty(SH),
    disp('compiling spont history of cell ...');
    clear SH;
    for ii=1:numel(DB),
        SH(ii) = structJoin(DB(ii), local_strip(JLgetRec(DB(ii).UniqueRecordingIndex)));
    end
    SH = sortAccord(SH, [SH.MinutesSinceRecStart]);
    putcache(CFN, 100, CP, SH);
end

if doPlot,
    set(figure,'units', 'normalized', 'position', [0.4 0.436 0.53 0.294]);
    S = JLdatastruct(Uidx);
    Tm = [SH.MinutesSinceRecStart];
    plot(Tm, [SH.SpontVar1], 'o-');
    xplot(Tm, [SH.FastSpontVar1], 'v-');
    xplot(Tm, [SH.SpontVar2], 'ro-');
    xplot(Tm, [SH.FastSpontVar2], 'rv-');
    xlabel('Rec time (minutes)');
    qq = [SH.SpontVar1];
    YL = max(ylim);
    ylim([0 median(qq)+min(YL, 4*iqr(qq))]);
    fenceplot(db.MinutesSinceRecStart, ylim, 'color', [0.4 1 0.4] ,'linewidth', 3);
    legend({'Spontvar 1' 'Spontvar 1 (detrended)' 'Spontvar 2' 'Spontvar 2 (detrended)' }, 'location', 'best');
    title(S.TTT);
end

% ==============
function db = local_strip(db);
[db.rec db.SPT db.SPTraw db.Tsnip db.Snip ] = deal([]);



