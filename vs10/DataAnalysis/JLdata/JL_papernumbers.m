function Y = JL_papernumbers(k)
% JL_papernumbers - metrics, stats, &c, for MSO paper
%    JL_papernumbers(k)
%    k=1: halfwidth of spont EPSPs
%    k=2: neurophonix
%      3: EPSP-AP latency stats
%      4: spont AP rates


global DB NNTP NNTPsh
if isempty(DB),
    DB = structMerge('UniqueRecordingIndex', JLdbase, JL_rec_select); 
    DB = DB([DB.CanBeUsed]);
end
if isempty(NNTP),
    NNTP = structMerge('UniqueRecordingIndex', JLdbase, JL_rec_select, JL_rec_select, JLcycleStats, JLcycleCorrs, JLanovaStats, JLspikeStats, ...
        retrieve_dbase('D:\processed_data\JL\NNTP\PhaseLock.dbase'));
    NNTP = NNTP([NNTP.CanBeUsed]); % restrict based on quality of rec and AP-detectability
    T_LAT = retrieve_dbase('D:\processed_data\JL\NNTP\JL_APana.dbase');
    NNTP = structMerge('UniqueRecordingIndex', NNTP, T_LAT);
end
NNTPsh = NNTP([NNTP.ShouldBeUsed]);
NNTPc = NNTP([NNTP.CanBeUsed]);
% unique([DB.icell_run])
%      1     4     6     7    13    14    15    16    17    18    19    20    22    23    25    30    31    32    33    35
%     37    38    40    41    42

switch k
    case 1, % halfwidth of spont EPSPs
        % manual
        %         JLspontPeaks3(Uidx,[],[],-1)
        %         DT = [];
        %         qq=ginput(2); DT = [DT, diff(qq(:,1))]; % repeat this line
        %         sprintf('%f  ', mean(DT), median(DT), std(DT), iqr(DT), numel(DT)) % for display
        %         save('D:\processed_data\JL\spontEPSP\DT_cellXXX', 'DT') % save DT
        % icell_run Uidx    MEAN     MEDIAN    STD       IQR      N   
        %  1   178013810 0.189992  0.188758  0.036356  0.041946  34
        %  4   191020904 0.262812  0.270321  0.046411  0.046607  36
        %  6   195034311 0.240292  0.243289  0.027324  0.029362  28
        %  7   196020512 0.183845  0.176174  0.027981  0.033977  35
        % 13   201010108 0.240419  0.239094  0.030918  0.041946  38
        % 14   201020404 0.239693  0.243289  0.026000  0.037752  35
        % 15   201031304 0.191781  0.186871  0.023427  0.027181  38
        % 16   204113005 0.255872  0.251678  0.029774  0.046607  40
        % 17   204126101 0.209994  0.211409  0.020154  0.026426  32
        % 18   204137001 0.248622  0.245386  0.035963  0.052852  14
        % 19   204147807 0.222370  0.224245  0.017692  0.022085  29
        % 20   209020401 0.208036  0.203859  0.030357  0.037752  47
        % 22   209041201 0.183725  0.181208  0.020195  0.030201  27
        % 23   209051703 0.203859  0.203859  0.030549  0.049077  28
        % 25   214010107 0.155677  0.158557  0.018115  0.026237  43
        % 30   214061603 0.203455  0.189550  0.039742  0.039490  44
        % 31   216010201 0.148269  0.149660  0.026775  0.044277  32
        % 32   216020801 0.140054  0.136326  0.019682  0.027265  36
        % 33   219010101 0.185582  0.190268  0.029234  0.035675  29
        % 37   219055105 0.252936  0.256711  0.029986  0.030201  34
        % 38   219065304 0.206068  0.201820  0.036250  0.039753  36
        % 40   223020404 0.177510  0.170630  0.031400  0.030273  24
        % 41   223030801 0.200976  0.203859  0.025085  0.035675  33
        % 42   223041205 0.208140  0.207936  0.035058  0.048926  30
        dd=dir('D:\processed_data\JL\spontEPSP\DT_cell*.mat');
        for ii=1:numel(dd),
            DDTT(ii) = load(fullfile('D:\processed_data\JL\spontEPSP\', dd(ii).name));
        end
        DT = [DDTT.DT];
        sprintf('%f  ', mean(DT), median(DT), std(DT), iqr(DT), numel(DT))
        %    MEAN     MEDIAN    STD       IQR      N   
        %  0.205981  0.203859  0.044485  0.061242  802
        Y.meanHalfwidth = mean(DT);
        Y.medianHalfwidth = median(DT);
        Y.stdHalfwidth = std(DT);
        Y.iqrHalfwidth = iqr(DT);
        Y.N = numel(DT);
        Y.Ncell = numel(unique([DB.icell_run]));
    case 2, % k=2: neurophonix
        DB0 = JLdbase;
        % cell 12 = 198/2: 2-18-BFS, 2-19-BFS, 2-20-BFS.
        % iseries_run:        38        40       43
        db12_18 = DB0([DB0.iseries_run]==38);
        db12_19 = DB0([DB0.iseries_run]==40);
        db12_20 = DB0([DB0.iseries_run]==43);
        nf12 = [db12_18 db12_19 db12_20];
        WS12 = JLcycleStats([nf12.UniqueRecordingIndex]);
        Y.meanNeuroPhonRMS = mean(sqrt([WS12.VarDriv]));
        Y.stdNeuroPhonRMS = std(sqrt([WS12.VarDriv]));
        Y.Nrec = numel(WS12);
    case 3, % EPSP-AP latency stats
        qq=retrieve_dbase('D:\processed_data\JL\NNTP\JL_APana.dbase');
        N = nansum([qq.Nlat]);
        meanLat = nansum([qq.MeanLat].*[qq.Nlat])/N;
        stdLat = sqrt(nansum([qq.StdLat].^2.*[qq.Nlat])/N);
        meanCVLat = nanmean([qq.StdLat]./[qq.MeanLat]);
        Y = CollectInStruct(N, meanLat, stdLat, meanCVLat);
    case 4, % spont AP rates
        AlliCell = unique(unique([NNTP.icell_run]));
        SpontRate = [];
        for ii=1:numel(AlliCell),
            nntp = NNTP([NNTP.icell_run]==AlliCell(ii));
            spr = mean([nntp.SpontRate1]);
            SpontRate = [SpontRate spr];
        end
        meanSpontRate = mean(SpontRate);
        medianSpontRate = median(SpontRate);
        stdSpontRate = std(SpontRate);
        iqrSpontRate = iqr(SpontRate);
        Y = CollectInStruct(SpontRate, meanSpontRate, stdSpontRate, medianSpontRate, iqrSpontRate);
    case 5, % variation of spont rate within series
        AllCell = unique([NNTPsh.icell]);
        for ii=1:numel(AllCell),
            nntp = NNTPsh([NNTPsh.icell]==AllCell(ii));
            Y.VarSpont1{ii} = [nntp.VarSpont1];
            Y.MeanSpontRMS(ii) = mean(sqrt([nntp.VarSpont1]));
            Y.StdSpontRMS(ii) = std(sqrt([nntp.VarSpont1]));
            Y.MinSpontRMS(ii) = min(sqrt([nntp.VarSpont1]));
            Y.MaxSpontRMS(ii) = max(sqrt([nntp.VarSpont1]));
        end
    case 6, % percentage of APs in analysis window
        DC = 1.5*[NNTPc.NspikeInStim]./[NNTPc.AnaDur];
        Y.median = median(100*DC);
        Y.iqr = iqr(100*DC);
    case 7, % EPSP-AP latencies
        LatDiff = [];
        LAT = retrieve_dbase('D:\processed_data\JL\NNTP\JL_APana.dbase'); % latency data
        C = retrieve_dbase('D:\processed_data\JL\JLcompareSpikeCounts\JLcompareSpikeCounts_10.dbase');  %triplets
        FNS = {'UniqueTripletIndex', 'ID___________________' 'iexp' 'icell' 'icell_run' 'freq' 'SPL'};
        for ii=1:numel(C),
            Lat_I = LAT([LAT.UniqueRecordingIndex]==C(ii).Uidx_I);
            Lat_C = LAT([LAT.UniqueRecordingIndex]==C(ii).Uidx_C);
            if isempty(Lat_I)
                warning(['Missing latency data' JLtitle(C(ii).Uidx_I)])
            elseif isempty(Lat_C)
                warning(['Missing latency data' JLtitle(C(ii).Uidx_C)])
            else,
                Uidx_I = Lat_I.UniqueRecordingIndex;
                Nlat_I = Lat_I.Nlat;
                MeanLat_I = Lat_I.MeanLat;
                StdLat_I = Lat_I.StdLat;
                Uidx_C = Lat_C.UniqueRecordingIndex;
                Nlat_C = Lat_C.Nlat;
                MeanLat_C = Lat_C.MeanLat;
                StdLat_C = Lat_C.StdLat;
                ld = CollectInStruct('-EPSP_AP_latency', Uidx_I, Nlat_I, MeanLat_I, StdLat_I ,Uidx_C, Nlat_C, MeanLat_C, StdLat_C);
                ld = structJoin(structpart(C(ii), FNS) ,ld);
                LatDiff = [LatDiff, ld];
            end
        end
        AllCells = unique([LatDiff.icell_run]);
        figure;
        set(gcf,'units', 'normalized', 'position', [0.47 0.619 0.248 0.28])
        Ncell = numel(AllCells);
        [meanlat_I, meanlat_C] = deal(nan(1,Ncell));
        for ii=1:Ncell,
            ld = LatDiff([LatDiff.icell_run]==AllCells(ii));
            meanlat_I(ii) = 1e3*nanmean([ld.MeanLat_I]);
            if ii==1,
                meanlat_I(1) = meanlat_I(1).*0.8; % correct miscounts in JL_APThr
            elseif ii==11,
                meanlat_I(11) = meanlat_I(11).*1.3; % correct miscounts in JL_APThr
            end
            meanlat_C(ii) = 1e3*nanmean([ld.MeanLat_C]);
            xplot(meanlat_I(ii), meanlat_C(ii), [ploco(ii) ploma(ii)], 'markersize', 8);
        end
        TT = ttest(denan(meanlat_I-meanlat_C))
        xlim([0 388]); ylim(xlim);
        xplot(xlim,xlim, 'k');
        Y = CollectInStruct(meanlat_I, meanlat_C);
        set(gca, 'fontsize', 12)
        xlabel('Ipsi EPSP-AP latency (\mus)', 'fontsize' ,13)
        ylabel('Contra EPSP-AP latency (\mus)', 'fontsize' ,13)
    otherwise, error(['Unkown NNTP k=' num2str(k)]);
end


