function Y = JL_paperfigs(Panel, varargin)
% JL_paperfigs - data selection for figures of MSO paper
%    1A: complex subthreshold events in beat responses
%        cell 30 (214/6) B 700 Hz, 50 dB. (irec = 157; Uidx=214062507)
%    1B: spontaneous subthreshold events
%    1C: subthr events and spikes in beat responses. V(t) and dV/dt
%        cell 30 (214/6) B 400 Hz, 50 dB. (irec = 154; Uidx=214062504)
%        Spike at t=739 ms
%    1D: histogram of downward slopes for spike detection
%    1E:  spike waveforms illustrating corr btw EPSP size & latency
%    1F:  EPSP size vs AP latency
%
%    2A:  tone-delay curve
%    2B:  averaging of stim cycles
%    2C:  2D plot of subthr inputs plus spikes
%    2D:  population data predicted vs measured best IPD
% 
%
%    3A:  cycle-avgs, freqs 200-1000, 80 dB.
%    3B:  freq dep bin tuning
%  3CDE:  3 bin receptive fields
%
%    4A:  lin pred of 2D versus 2D data
%    4B:  population data: hist of VarAcc linear pred cf 3A
%    4C:  cycle-averaged waveforms mon / bin & ipsi /contra
%    4D:  2D waveforms predicted from monaural responses
%    4E:  2D variance
%
%    5A:  waveforms I,C,B showing thresholding
%    5B:  hist of ratio binaural/(summed monaural) spike rates
%    5C:  power-law example
%
%    S2F: ITD population data
%    S3A: time spent below spont var as a fnc of freq
%
%    ====SOM figures=====
%
%    102: best ITD varies w freq
%
%    PRINT
%     print(gcf, '-r400', '-dpng', '\\Ee1285a\MSO$\Figure_X.png')

if nargin>1,
    JL_paperfigs(Panel);
    for ii=1:nargin-1,
        JL_paperfigs(varargin{ii});
    end
    return;
end

global Usbu % idenitifiers of recordings tagged "should be used"
if isempty(Usbu),
    SL = JL_rec_select;
    Usbu = [SL([SL.ShouldBeUsed]).UniqueRecordingIndex];
end

FZ = 'fontsize';
MFC = 'markerfacecolor';
TraceLW = 2.5; % line width for raw traces
switch upper(Panel)
    case 'PRINTALL'
        JL_paperfigs 1
        print(gcf, '-r800', '-dpng', '\\Ee1285a\MSO$\Figure_1.png');
        JL_paperfigs 2
        print(gcf, '-r600', '-dpng', '\\Ee1285a\MSO$\Figure_2.png');
        JL_paperfigs 3
        print(gcf, '-r600', '-dpng', '\\Ee1285a\MSO$\Figure_3.png');
        JL_paperfigs 4
        print(gcf, '-r600', '-dpng', '\\Ee1285a\MSO$\Figure_4.png');
        pause(2);
        aa;
        for ifig = 102:111,
            str = num2str(ifig);
            JL_paperfigs(str);
            print(gcf, '-r400', '-dpng', ['\\Ee1285a\MSO$\Figure_' str ' .png']);
        end
    case '1' %
        h = local_fig(1);
        if ~isempty(h),delete(h); end
        JL_paperfigs 1A 1C 1D 1E 1F
        saveas(gcf, '\\Ee1285a\MSO$\Figure_1.fig');
    case '2'
        h = local_fig(2);
        if ~isempty(h),delete(h); end
        JL_paperfigs 2A 2B 2C 2D
        saveas(gcf, '\\Ee1285a\MSO$\Figure_2.fig');
    case '3'
        h = local_fig(2);
        if ~isempty(h),delete(h); end
        JL_paperfigs 3A 3B 3CDE 3F
        saveas(gcf, '\\Ee1285a\MSO$\Figure_3.fig');
    case '4'
        h = local_fig(3);
        if ~isempty(h),delete(h); end
        JL_paperfigs 4A 4B 4C 4D 4E 4F
        saveas(gcf, '\\Ee1285a\MSO$\Figure_5.fig');
    case '5'
        h = local_fig(4);
        if ~isempty(h),delete(h); end
        JL_paperfigs 5A 5B 5C
        saveas(gcf, '\\Ee1285a\MSO$\Figure_5.fig');
    case '1A' % complex events in beat responses
        Uidx = 214062507; % Exp 214 cell 6  --- B  700 Hz   50 dB  (series 7)
        local_fig(1);
        Xmin = 0.07; Y0 = 0.40;
        ha1 = axes('pos', [Xmin Y0 0.35 0.52]);
        [W, Rc] = JLgetRec(Uidx);
        dplot(W.dt, Rc(1:740000), 'k', 'linewidth', 1);
        xplot([500 6500], [1 1], 'color', [0.2 0.2 1], 'linewidth', 5)
        ylim([-10 7]);
        xlim([-400 9500]);
        [hb, htx, hty] = local_scalehook([8300 2.5], 3, -1000, '1 s', 0.2, 1, '1 mV', 100, 12);
        xlabel('Time (ms)');
        ylabel('V_{rec} (mV)')
        local_mark(ha1, 'A');
        % zoom in to see subthr events and spikes in beat responses. V(t) and dV/dt
        Twin = [4302, 4312];
        xplot([Twin(1)*[1 1] 144*[1 1]], [4 2 -1.3 -1.8], '--', 'color', 0.6*[1 1 1], 'linewidth', 2)
        xplot([Twin(2)*[1 1] 8000*[1 1]], [4 2 -1.3 -1.8], '--', 'color', 0.6*[1 1 1], 'linewidth', 2)
        ha2 = axes('pos', [Xmin+0.02 Y0+0.1428 0.33 0.15]);
        [W, Rc] = JLgetRec(Uidx);
        Tm = Xaxis(Rc, W.dt);
        Rc = Rc(betwixt(Tm, Twin));
        Tm = Tm(betwixt(Tm, Twin));
        plot(Tm, Rc, 'k', 'linewidth', TraceLW);
        xlim(Twin+[0 2]);
        ylim([2.1 5.2]);
        dotColor = [1 0.1 0.1];
        xplot([4304.3 4310.04], [4.7 4.9], 'o', 'markerfacecolor', dotColor, 'color', dotColor, 'markersize' ,7)
        xlabel('Time (ms)');
        ylabel('dV/dt (V/s)      V_{rec} (mV)')
        [hb, htx, hty] = local_scalehook([4312.5 3], 3, -0.001, '', 0.1, 1, '1 mV', 0.2, 12);
        local_mark(ha2, 'B', -0.04, -0.07);
        % time derivative
        ha3 = axes('pos', [Xmin+0.02 Y0 0.33 0.15]);
        [W, Rc] = JLgetRec(214062507);
        Tm = Xaxis(Rc, W.dt);
        Twin = [4302-10, 4312+10];
        Rc = Rc(betwixt(Tm, Twin));
        Tm = Tm(betwixt(Tm, Twin));
        dRc = diff(smoothen(Rc, 0.08, W.dt))/W.dt;
        Twin = [4302, 4312];
        dRc = dRc(betwixt(Tm, Twin));
        Tm = Tm(betwixt(Tm, Twin));
        ScaleFac = 1/8;
        xplot(Tm, 1.5+ScaleFac*dRc, 'k', 'linewidth', TraceLW)
        [hb, htx, hty] = local_scalehook([4312.5 0.7], 3, -1, '1 ms', 0.1, 10*ScaleFac, '10 V/s', 0.2, 12);
        xlim(Twin+[0 2]);
        ylim([-0.5 3]);
        %local_mark(h2, 'C');
        set([ha1 ha2 ha3], 'vis', 'off');
    case '1C' % spontaneous subthreshold events
        local_fig(1);
        ha = axes('pos', [0.09  0.20  0.33  0.15]);
        [W, Rc] = JLgetRec(214062507);
        Tm = Xaxis(Rc, W.dt);
        Twin = [1120, 1130];
        Rc = Rc(betwixt(Tm, Twin));
        Tm = Tm(betwixt(Tm, Twin));
        plot(Tm, Rc, 'k', 'linewidth', TraceLW);
        [hb, htx, hty] = local_scalehook([Twin(2)+0.5 3], 3, -1, '1 ms', 0.1, 0.93, '1 mV', 0.2, 12);
        %[hb, htx, hty] = local_scalehook([Twin(2)+1.3 2.93], 3, -0.001, '', 0.1, 1, '1 mV', 0.2, 12);
        xlim(Twin+[0 2]);
        ylim(2+[0 3.1]); % has to match ylim span in ha2 of panel A
        xlabel('Time (ms)');
        ylabel('V_{rec} (mV)')
        set(ha,'visib', 'off');
        local_mark(ha, 'C');
    case '1D' % histogram of downward slopes for spike detection
        local_fig(1);
        ha = axes('pos', [0.55 0.75  0.25  0.17]);
        TH = JLspikeThr(214062506, nan, 0);
        local_APthrplot(TH);
        set(ha,'box', 'off')
        local_mark(ha, 'D');
    case '1E' % spike waveforms illustrating corr btw EPSP size & latency
        local_fig(1);
        ha = axes('pos', [0.55 0.375 0.30 0.30]);
        %ha = axes('pos', [0.1 0.8 0.35 0.15]);
        Uidx = 214062204;
        A=JL_APana(Uidx,0); % 0=no plot
        W = JLwaveforms(Uidx);
        SNP = W.Snip(:,W.APinStim);
        [Lat, iord] = sort(A.Lat);
        SNP = SNP(:,iord);
        Vepsp = A.Vepsp(iord);
        %                         plot(Vepsp, Lat, '.w')
        %                         for ii=1:numel(Vepsp),
        %                             text(Vepsp(ii), Lat(ii), num2str(ii))
        %                         end
        %                         pause; liuhljhlgbvds
        iord(48)
        igrow = fliplr([52 45 40 29 14 5 3]);
        SNP = SNP(:, igrow);
        Tsn = bsxfun(@plus, W.Tsnip, Lat(igrow));
        Tsn(:,1) = Tsn(:,1) + 0.08;
        Tsn(:,2) = Tsn(:,2) + 0.10;
        xplot(Tsn+pmask(Tsn<=1.2), SNP, 'linewidth',2);
        xlim([-0.36 1.5]);
        ylim([0.5 5.3]);
        %scalebar
        %[hb, ht] = local_scalebar(0, 2, 0.2, 0.1, 2 ,'0.2 ms', 0.2, 11);
        [hb, htx, hty] = local_scalehook([1.3 1.6], 3, -0.2, '0.2 ms', 0.2, 1, '1 mV', 0.05, 12);
        text(0, 2.8, '\uparrow', 'fontsize', 20 ,'fontweight', 'bold', 'horizontalalign', 'center')
        text(0, 2.2, 'EPSP', 'fontsize', 13 ,'fontweight', 'normal', 'horizontalalign', 'center')
        %text(0.36, 4.7, '\leftarrow', 'fontsize', 20 ,'fontweight', 'bold', 'horizontalalign');
        Xar = [0.2358  0.5206  0.37]-0.05; Yar = [ 4.4989 4.7239  3.5]+0.1;
        xplot(Xar, Yar, 'k', 'linewidth', 2);
        xplot(Xar(1)+[0.05 0 0.06], Yar(1)+[0.1 0 0], 'k', 'linewidth', 2);
        xplot(Xar(3)+[0.015 0 0.04], Yar(3)+[0.26 0 0.2], 'k', 'linewidth', 2);
        text(Xar(2), Yar(2), 'APs', 'fontsize', 13 ,'fontweight', 'normal', 'backgroundcol', 'w');
        set(gca,'visib', 'off');
        local_mark(ha, 'E');
    case '1F' % EPSP size vs AP latency
        local_fig(1);
        ha = axes('pos', [0.55    0.16    0.3000    0.2406]);
        Uidx = 214062200;
        for ii=1:9,
            A(ii)=JL_APana(Uidx+ii,0); % 0=no plot
        end
        Ve = [A.Vepsp]; Lat = [A.Lat];
        iok = ~isnan(Lat);
        [Ve, Lat] = deal(Ve(iok), Lat(iok));
        [Ve, Lat] = sortAccord(Ve, Lat, Ve);
        dsize(Ve, Lat);
        P = polyfit(Ve(:), Lat(:),1);
        plot(Ve-3 ,Lat, '.', 'markersize', 14);
        %xplot(Ve-3, polyval(P, Ve), 'r', 'linewidth',2)
        xlim([0 1.9])
        set(ha,'box', 'off', 'fontsize', 12);
        xlabel('EPSP amplitude (mV)', 'fontsize', 14);
        ylabel('EPSP-AP delay (ms)', 'fontsize', 14);
        local_mark(ha, 'F');
    case '2A', % tone-delay curve
        local_fig(2);
        ha2 = axes('position', [0.1300    0.7221    0.15    0.2029], ...
            'YAxisLocation', 'right', 'xtick', []);
        ha = axes('position', [0.1300    0.7221    0.15    0.2029]);
        Uidx = 214062507; % Exp 214 cell 6  --- B  700 Hz   50 dB  (series 7)
        D = JLdbase(Uidx); % general info on the recording & stimulus
        S = JLdatastruct(Uidx); % specific info on stimulus etc
        W = JLwaveforms(Uidx); % waveforms and duration params
        T = JLcycleStats(Uidx); % cycle-based means & variances
        A = JL_APana(Uidx,0); % EPSP-AP latencies
        AD = JLanovaDetails(Uidx);
        TT = JLtitle(Uidx);
        dt = W.dt;
        SPT = W.SPTraw(:, W.APinStim)-T.t_stimOffset;
        Nspike = numel(SPT);
        Tbeat = 1e3/S.Fbeat; % beat period in ms
        T_I = 1e3/S.Freq1; % ipsi cycle in ms
        T_C = 1e3/S.Freq2; % contra cycle in ms
        SPT_I = mod(SPT, T_I)/T_I;
        SPT_C = mod(SPT, T_C)/T_C;
        SPT_B = mod(SPT, Tbeat)/Tbeat;
        Nhist = 18;
        Ph = ((1:Nhist)-0.5)/Nhist;
        [Nsp] = hist(SPT_B, Ph);
        Ph = mod(Ph+0.5,1)-0.5;
        ITD = T_I*Ph;
        [ITD, Nsp] = sortAccord(ITD, Nsp, ITD);
        AC = @(c)struct('color',c, 'markerfacecolor',c);
        % subthr
        axes(ha2);
        S = JLpeakPhase(Uidx);
        SAP = S.ShiftAddPeak;
        Nsap = numel(SAP);
        SAP = circshift(SAP,[round(Nsap/2) 1]).';
        SAP = SAP-min(SAP(:));
        ITDs = 1e3*linspace(-T_I/2, T_I/2, Nsap);
        xplot(ITDs, SAP, 'linewidth', 3, 'color', [0 0.8 0]);
        xlim(1e3*T_I*[-1 1]/2);
        ylim([0 0.599]);
        hyl = ylabel('Peak potential (mV)', FZ, 13, 'color', [0 0.8 0], 'rotation', -90);
        PYL = get(hyl, 'posi')
        set(hyl, 'posi', PYL+[170 0 0]);
        set(ha2,'ytick', 0:0.2:0.8);
        % spikes
        axes(ha);
        plot(ITD*1e3, Nsp, '-o', 'linewidth',2,'markersize', 7, AC('k'));
        xlim(1e3*T_I*[-1 1]/2);
        ylim([-2.5 28]);
        set(ha,'box', 'off', 'fontsize', 13, 'linewidth',2, 'color', 'none');
        set(ha2, 'box', 'off', 'fontsize', 13, 'linewidth',2);
        xlabel('ITD (\mus)', 'fontsize', 13);
        ylabel('Spike count', 'fontsize', 13);
        text(210, 24, '\downarrow', 'fontsize', 16, 'fontweight', 'bold', 'rotation', -20);
        local_mark(ha, 'A', 0.04);
    case '2B', % averaging of stim cycles
        local_fig(2);
        Xleft = 0.43;
        Ybt = 0.69;
        Uidx = 214062507; % Exp 214 cell 6  --- B  700 Hz   50 dB  (series 7)
        D = JLdbase(Uidx); % general info on the recording & stimulus
        S = JLdatastruct(Uidx); % specific info on stimulus etc
        W = JLwaveforms(Uidx); % waveforms and duration params
        T = JLcycleStats(Uidx); % cycle-based means & variances
        A = JL_APana(Uidx,0); % EPSP-AP latencies
        AD = JLanovaDetails(Uidx);
        TT = JLtitle(Uidx);
        dt = W.dt;
        Tipsi = 1e3/S.Freq1;  Tcontra = 1e3/S.Freq2;
        M = 60; DY = 0.6; Nsel = 15; tshift=0.7+[0 250];
        [mWi, selWi, dti] = local_cyclesnips(W, M, Tipsi, Nsel, DY, tshift(1));
        [mWc, selWc, dtc] = local_cyclesnips(W, M, Tipsi, Nsel, DY, tshift(2));
        ha1 = axes('position', [Xleft    Ybt    0.083    0.23]);
        local_stackplot(dti, mWi, selWi, struct('color', [0 0 0.8]), 'IPSI: 700 Hz');
        ha2 = axes('position', [Xleft+0.09 Ybt  0.083  0.23]);
        local_stackplot(dtc, mWc, selWc, struct('color', [0.8 0 0]), 'CONTRA: 704 Hz');
        local_mark(ha1, 'B', 0.04, 0.17, 'backgroundcolor', 'none');
        %         dplot(dti, selWi, 'k', 'linewidth',2);
        %         xdplot(dti, -1+0*mWi, 'k:', 'linewidth', 3);
        %         xdplot(dti, -4+4*mWi, 'b', 'linewidth', 3);
        %         xdplot(dti, 0.74*Nsel+Stim, 'b', 'linewidth', 3);
        %         text(0.2, 0.87*Nsel, 'IPSI: 700 Hz', 'color', 'b')
        %         xlim([0 Tipsi]);
        %         ylim([-5 0.97*Nsel]);
    case '2C', % 2D plot of subthr inputs plus spikes
        local_fig(2);
        Uidx = 214062504; % Exp 214 cell 6  --- B  400 Hz   50 dB  (series 7)
        T = JLanova2(Uidx);
        S = JLdatastruct(Uidx);
        ha = axes('position', [0.1300    0.3    1.3*[0.15    0.2029]], FZ, 12);
        Phase = linspace(0,1-1/T.NsamPerCycle,T.NsamPerCycle);
        %         imaxI = local_findpeak(T.ipsimean);
        %         imaxC = local_findpeak(T.contramean);
        %         nshift = round(T.NsamPerCycle/2)-[imaxI,imaxC];
        nshift = [0 0];%*nshift;
        mean2d = circshift(T.mean2d, nshift);
        linpred = circshift(T.linpred, nshift);
        var2d = circshift(T.var2d, nshift);
        Nspike = circshift(T.Nspike2d, nshift);
        Spike_phase_C = rem(T.Spike_phase_C+2+nshift(1)/T.NsamPerCycle, 1);
        Spike_phase_I = rem(T.Spike_phase_I+2+nshift(2)/T.NsamPerCycle, 1);
        %==
        CCC = contourf(Phase, Phase, mean2d);
        axis square
        xlabel('Contra phase (cycle)', FZ, 13)
        ylabel('Ipsi phase (cycle)', FZ, 13)
        dddx = 0.0085; dddy = dddx*0.7;
        xplot(Spike_phase_C+dddx*diff(xlim), Spike_phase_I+dddy*diff(ylim), '.k', 'markersize',8);
        xplot(Spike_phase_C, Spike_phase_I, '.w', 'markersize',8);
        set(ha,  FZ, 12);
        local_mark(ha, 'C', -0.1, 0.2, 'backgroundcolor', 'none');
    case '2D' % predicted versus measured best IPDs
        local_fig(2);
        ha = axes('position', [0.45  0.32    0.97*[0.15    0.22]]);
        qqq = local_IPDstuff;
        h = plot([qqq.predBestIPD], [qqq.measBestIPD], 'k.');
        IDpoints(h, '', [qqq.UniqueRecordingIndex], @(dum,Uidx)[num2str(Uidx) ' ' num2str(JL_icell_run(Uidx))], ...
            '2D', @(dum,Uidx)JLanova2(Uidx), ...
            'Cycles', @(dum,Uidx)JLbeatVar2(Uidx), ...
            'view', @(dum,Uidx)JLviewrec(Uidx), ...
            'spikes', @(dum,Uidx)JLspikes(Uidx)            );
        R = corrcoef([qqq.predBestIPD], [qqq.measBestIPD]);
        R = R(1,2);
        N = numel(qqq);
        xlim([-0.5 0.5]);
        ylim(xlim)
        xplot(xlim, xlim, 'k-');
        set(ha,  FZ, 13, 'box', 'off');
        xlabel('Predicted best IPD (cycle)', FZ, 13);
        ylabel('Measured best IPD (cycle)', FZ, 13);
        text(-0.4, 0.46, sprintf('{\\it R} = %1.2f \n({\\itN} = %d)', R, N), FZ, 12);
        Y = qqq;
        local_mark(ha, 'D', -0.2, 0.3, 'backgroundcolor', 'none');
    case '3A', % cycle-avgs, freqs 200-1000, 80 dB.
        local_fig(3);
        Xleft = 0.12;
        Ybt = 0.63;
        ha = axes('position', [Xleft-0.0001    Ybt    0.12    0.28]);
        ha2 = axes('position', [Xleft    Ybt    0.12    0.28], 'YAxisLocation', 'right');
        ioff = [2 3 4 6 8 10]
        Uidx = 214062200+ioff;
        freq = 100*ioff;
        NN = numel(Uidx);
        DY = 0.75;
        Cdelay = 5.5; % ms delay
        for ii=1:NN,
            W = JLwaveforms(Uidx(ii));
            [wi, wc] = deal(W.IpsiMeanrec, W.ContraMeanrec);
            [wi, wc] = deal(wi-mean(wi(:)), wc-mean(wc(:)));
            [Ni, Nc] = deal(numel(wi),numel(wc));
            fr = freq(ii)/1e3; Dphi = Cdelay*fr; Nshift = round((0.5-Dphi)*[Ni Nc]);
            [wi, wc] = deal(circshift(wi,Nshift(1)), circshift(wc,Nshift(2)));
            xdplot([1/Ni -0.5], (ii-1)*DY+wi, ipsicolor, 'linewidth', 2);
            xdplot([1/Nc -0.5], (ii-1)*DY+wc, contracolor, 'linewidth', 2);
        end
        axes(ha);
        set(ha, 'ytick', DY*(0:NN-1), 'yticklabel', Words2cell(num2str(freq)));
        YL = DY*[-0.65 NN+0.2];
        ylim(YL)
        ylabel('Ipsi Frequency (Hz)', 'fontweight', 'normal', 'color', getfield(ipsicolor,'color'), FZ,13);
        xlabel('Phase (cycle)', 'fontweight', 'normal', FZ, 13);
        %hleg = legend({'Ipsi' 'Contra'})
        set(ha2, 'ytick', DY*(0:NN-1), 'yticklabel', Words2cell(num2str(freq+4)), 'xtick',[]);
        axes(ha2);
        ylim(YL)
        hyl = ylabel('Contra Frequency (Hz)', 'fontweight', 'normal', 'color', ...
            getfield(contracolor,'color'), FZ,13, 'rotation' ,-90);
        PYL = get(hyl, 'posi')
        set(hyl, 'posi', PYL+[0.23 0 0]);
        set([ha ha2],  FZ, 13);
        %local_scalebarY(X,Y, DX, DY, Lw ,Txt, dxt, Fsz);
        local_scalebarY(-0.14, (NN-0.38)*DY, 0.05, 0.5, 2 ,'0.5 mV', 0.05, 10)
        local_mark(ha2, 'A', -0.3, 0.17, 'backgroundcolor', 'none');
    case '3B' % freq dep bin tuning
        local_fig(3);
        ha = axes('position', [0.42  0.58    0.15    0.37]);
        S = JLpeakPhase(-214062507);
        SAP = [S.ShiftAddPeak].';
        Freq = [S.Freq_ipsi].';
        Nph = size(SAP,2);
        Nsh = round(Nph/2);
        SAP = circshift(SAP, [0 Nsh]);
        PH = linspace(-0.5,0.5,Nph);
        dsize(PH, Freq, SAP)
        %contourf(PH, Freq, SAP);
        S = S([2 3 4 6 8 10]);
        set(ha,  FZ, 12);
        DY = 0.6;
        NN = numel(S);
        AC = @(c)struct('color',c, 'markerfacecolor',c);
        for ii=1:NN,
            % subthr
            sap = S(ii).ShiftAddPeak.';
            sap = circshift(sap, [0 Nsh]);
            sap = sap-min(sap);
            maxsap = max(sap);
            Yoffset = (ii-1)*DY;
            xplot(PH, Yoffset + sap, 'linewidth', 3, 'color', [0 0.8 0]);
            xplot(PH, Yoffset + 0*sap, 'linewidth', 1, 'color', [0 0 0]);
            % spikes
            Uidx = S(ii).UniqueRecordingIndex;
            D = JLdbase(Uidx); % general info on the recording & stimulus
            DS = JLdatastruct(Uidx); % specific info on stimulus etc
            W = JLwaveforms(Uidx); % waveforms and duration params
            T = JLcycleStats(Uidx); % cycle-based means & variances
            A = JL_APana(Uidx,0); % EPSP-AP latencies
            AD = JLanovaDetails(Uidx);
            TT = JLtitle(Uidx);
            dt = W.dt;
            SPT = W.SPTraw(:, W.APinStim)-T.t_stimOffset;
            Nspike = numel(SPT);
            Tbeat = 1e3/DS.Fbeat; % beat period in ms
            T_I = 1e3/DS.Freq1; % ipsi cycle in ms
            T_C = 1e3/DS.Freq2; % contra cycle in ms
            SPT_I = mod(SPT, T_I)/T_I;
            SPT_C = mod(SPT, T_C)/T_C;
            SPT_B = mod(SPT, Tbeat)/Tbeat;
            Nhist = 12;
            Ph = ((1:Nhist)-0.5)/Nhist;
            [Nsp] = hist(SPT_B, Ph);
            Ph = mod(Ph+0.5,1)-0.5;
            ITD = T_I*Ph;
            [Ph, Nsp] = sortAccord(Ph, Nsp, Ph);
            Nsp = maxsap*Nsp/max(Nsp);
            %
            Yscale = 0.05;
            xplot(Ph, Yoffset + Nsp, '-o', AC('k'), 'markersize',5, 'linewidth',2);
            %xplot(SPH, dfreq+Freq(ii), '.w', 'markersize',5);
        end
        %==
        [hb, ht] = local_scalebarY(-0.44, (NN-0.38)*DY, 0.05, 0.5, 2 ,'0.5 mV', 0.05, 12);
        set(hb, 'color', [0 0.8 0]);
        set(ht, 'color', [0 0.8 0]);
        ylim([-0.2 (NN+0.2)*DY]);
        freq = [S.Freq_ipsi];
        set(gca, 'ytick', DY*(0:NN-1), 'yticklabel', Words2cell(num2str(freq)), FZ, 13);
        xlabel('IPD (cycle)', FZ, 13);
        ylabel('Ipsi fqrequency (Hz)', FZ, 13);
        local_mark(ha, 'B', -0.3429, 0.0601, 'backgroundcolor', 'none');
    case '3CDE' % 1/3 more bin receptive fields 
        local_fig(3);
        LM = @(ha,CH)local_mark(ha, CH,-0.05,0.2);
        ha = axes('position', [0.68  0.58    0.15    0.34]);
        Uidx = 204112801; JLtitle(Uidx)
        [S, mWvI, mWvC] = JLpeakPhase(-Uidx, 5.5);
        local_BRF_plot(S, mWvI, mWvC);
        LM(gca, 'C')
        ha = axes('position', [0.1  0.1   0.15    0.34]);
        Uidx = 214061804; JLtitle(Uidx)
        [S, mWvI, mWvC] = JLpeakPhase(-Uidx, 5.5);
        local_BRF_plot(S, mWvI, mWvC,0);
        LM(gca, 'D')
        ha = axes('position', [0.4  0.1   0.15    0.34]);
        Uidx = 219065301; JLtitle(Uidx)
        [S, mWvI, mWvC] = JLpeakPhase(-Uidx, 5.5);
        local_BRF_plot(S, mWvI, mWvC,0);
        LM(gca, 'E')
        %=
    case '3F' % theoriticial bin rec fld
        local_fig(3);
        ha = axes('position', [0.68  0.1    0.15    0.34]);
        LM = @(ha,CH)local_mark(ha, CH,-0.05,0.2);
        Uidx = 219065301; JLtitle(Uidx)
        [S, mWvI, mWvC] = JLpeakPhase(-Uidx, 5.5);
        %=
        Nph = size(mWvI,1); % # phase values
        Nh = round(Nph/2);
        dph = 1/Nph;
        PH = dph*(0:Nph-1).'-0.5;
        [Wmin, Wmax] = minmax([mWvI(:); mWvC(:)]);
        maxFreq = max([S.Freq_ipsi]);
        ww2 = 0.4^2;
        for ifr=1:numel(S),
            freq = S(ifr).Freq_ipsi;
            Wfr = exp(-((freq-500)/500).^2);
            best_ph = 0.3*freq/maxFreq;
            fxp = @(PH)exp(-(PH-best_ph).^2/ww2);
            xp = fxp(PH)+fxp(PH-1)+fxp(PH+1);
            xp = xp-min(xp);
            xp = Wfr*xp/max(xp);
            sap = Wmin + 0.9*(Wmax-Wmin)*xp;
            sap = circshift(sap ,Nh);
            S(ifr).ShiftAddPeak = sap;
        end
        local_BRF_plot(S, mWvI, mWvC,0);
        LM(gca, 'F')
    case '4A', % lin pred of 2D versus 2D data
        local_fig(3);
        ha1 = axes('position', [0.1  0.6    1.1*[0.15    0.22]], FZ, 12);
        ha2 = axes('position', [0.28  0.6    1.1*[0.15    0.22]], FZ, 12);
        %Uidx = 214062202;
        Uidx = 214061903;
        T = JLanova2(Uidx);
        S = JLdatastruct(Uidx);
        Phase = linspace(0,1-1/T.NsamPerCycle,T.NsamPerCycle);
        %imaxI = local_findpeak(T.ipsimean);
        %imaxC = local_findpeak(T.contramean);
        %nshift = round(T.NsamPerCycle/2)-[imaxI,imaxC];
        nshift = [0 0];%*nshift;
        %=measured
        mean2d = circshift(T.mean2d, nshift);
        linpred = circshift(T.linpred, nshift);
        Nspike = circshift(T.Nspike2d, nshift);
        mean2d = mean2d-mean(mean2d(:));
        linpred = linpred-mean(linpred(:));
        VarAcc_LIN = VarAccount(mean2d, linpred)
        Spike_phase_C = rem(T.Spike_phase_C+2+nshift(1)/T.NsamPerCycle, 1);
        Spike_phase_I = rem(T.Spike_phase_I+2+nshift(2)/T.NsamPerCycle, 1);
        axes(ha1);
        CCC = contourf(Phase, Phase, mean2d);
        axis square
        hx = xlabel('Contra phase (cycle)', FZ, 13);
        set(hx, 'position', [1.05   -0.1560    1.0001]);
        ylabel('Ipsi phase (cycle)', FZ, 13)
        SpikeMarkSize = 10;
        dddx = 0.0085; dddy = dddx*0.7;
        xplot(Spike_phase_C+dddx*diff(xlim), Spike_phase_I+dddy*diff(ylim), '.k', 'markersize',SpikeMarkSize);
        xplot(Spike_phase_C, Spike_phase_I, '.w', 'markersize',SpikeMarkSize);
        title('Recording', 'fontangle', 'italic')
        ca = caxis;
        %=predicted
        axes(ha2);
        contourf(Phase, Phase, linpred);
        caxis(ca);
        axis square
        set(ha2, 'yticklabel','');
        title('Linear prediction', 'fontangle', 'italic')
        local_mark(ha1, 'A', -0.2, 0.3, 'backgroundcolor', 'none');
    case '4B' % population data: hist of VarAcc linear pred cf 3A
        local_fig(3);
        ha = axes('position', [0.53  0.6    0.1   0.24], FZ, 12);
        DBFN = 'D:\processed_data\JL\JLanova2\AnovaStats.dbase';
        Tst = retrieve_dbase(DBFN);
        qqq = local_IPDstuff();
        Tst = Tst(ismember([Tst.UniqueRecordingIndex], [qqq.UniqueRecordingIndex]));
        UU = [Tst.UniqueRecordingIndex]; ddd=JLdbase(UU);
        Ncell = numel(unique([ddd.icell_run]))
        hist([Tst.Anova_VarAcc], linspace(0,100,200));
        xlim([75, 100]);
        ylim([0 125]);
        xlabel('Var accounted for (%)', FZ, 13);
        ylabel('# Recordings per bin', FZ, 13);
        text(83, 110, ['{\itN}=' num2str(numel(Tst))], FZ, 12);
        local_mark(ha, 'B', -0.4, 0.3, 'backgroundcolor', 'none');
    case '4C' % cycle-averaged waveforms mon / bin & ipsi /contra
        local_fig(3);
        ha = axes('position', [0.7  0.6    1.1*[0.15    0.22]], FZ, 12);
        %Uidx_B = 214062202;
        %Uidx_I = 214062402;
        %Uidx_C = 214062302;
        Uidx_B = 214061903;
        Uidx_I = 214062103
        Uidx_C = 214062003
        W = JLwaveforms(Uidx_B);
        Wbi = W.IpsiMeanrec;
        Wbc = W.ContraMeanrec;
        W = JLwaveforms(Uidx_I);
        Wmi = W.IpsiMeanrec;
        W = JLwaveforms(Uidx_C);
        Wmc = W.ContraMeanrec;
        II = struct('linewidth', 2, 'color', [0 0 0.85]);
        CC = struct('linewidth', 2, 'color', [0.85 0 0]);
        DY = 1.4;
        local_cplot(0, DY, Wbi, II);
        local_cplot(1.1, DY, Wbc, CC);
        local_cplot(0, 0, Wmi, II);
        local_cplot(1.1, 0 ,Wmc, CC);
        xlim([-0.05 2.15]);
        ylim([-0.6 2.15]);
        set(ha, 'box', 'on', 'xtick', [], 'ytick', [])
        xplot(xlim, mean(ylim)+0*ylim, 'k');
        xplot(mean(xlim)+0*xlim, ylim, 'k');
        text(0.25, -0.1, 'IPSI', 'units', 'normalized', 'horizontalalign', 'center', FZ, 13, 'color', II.color);
        text(0.75, -0.1, 'CONTRA', 'units', 'normalized', 'horizontalalign', 'center', FZ, 13, 'color', CC.color);
        text(-0.15, 0.25, 'MON', 'units', 'normalized', 'horizontalalign', 'center', FZ, 13, 'color', 'k', 'rotation', 0);
        text(-0.15, 0.75, 'BIN', 'units', 'normalized', 'horizontalalign', 'center', FZ, 13, 'color', 'k', 'rotation', 0);
        %                       local_scalehook(XY, Lw, DX, tX, dxt,   DY, tY, dyt, FZ);
        %local_scalebarY(0.2,-0.4, DX, DY, Lw ,Txt, dxt, Fsz);
        local_scalebarY(0.14,-0.36, 0.1, 0.3, 2 ,'0.3 mV', 0.05, 10);
        local_mark(ha, 'C', -0.2, 0.3, 'backgroundcolor', 'none');
    case '4D' % 2D waveforms predicted from monaural responses
        local_fig(3);
        ha = axes('position', [0.1  0.2    1.1*[0.15    0.22]], FZ, 12);
        %Uidx_B = 214062202;
        %Uidx_I = 214062402;
        %Uidx_C = 214062302;
        Uidx_B = 214061903;
        Uidx_I = 214062103
        Uidx_C = 214062003
        T = JLanova2(Uidx_B);
        nshift = [0 0];%*nshift;
        BIN = circshift(T.mean2d, nshift);
        BIN = BIN-mean(BIN(:));
        T = JLanova2(Uidx_I);
        PHi = linspace(0,1-1/T.NsamPerCycle,T.NsamPerCycle);
        Wmi = circshift(T.mean2d, nshift);
        Wmi = mean(Wmi,2);
        Wmi = Wmi-mean(Wmi(:));
        T = JLanova2(Uidx_C);
        PHc = linspace(0,1-1/T.NsamPerCycle,T.NsamPerCycle);
        nshift = [0 0];%*nshift;
        Wmc = circshift(T.mean2d, nshift);
        Wmc = mean(Wmc,1);
        Wmc = Wmc-mean(Wmc(:));
        Wmc = Wmc-mean(Wmc(:));
        Ni = numel(Wmi); PHi = (0:Ni-1)/Ni;
        Nc = numel(Wmc); PHc = (0:Nc-1)/Nc;
        [Wmc, Wmi] = SameSize(Wmc(:).', Wmi(:));
        linpred = Wmi+Wmc;
        dsize(PHi, PHc, linpred)
        contourf(PHi, PHc, linpred);
        VARAK = VarAccount(BIN(:), linpred(:));
        axis square
        title('Lin. pred. (monaural)', 'fontangle', 'italic')
        xlabel('Contra phase (cycle)', FZ, 13);
        ylabel('Ipsi phase (cycle)', FZ, 13)
        local_mark(ha, 'D', -0.2, 0.3, 'backgroundcolor', 'none');
    case '4E', % 2D variance
        local_fig(3);
        ha = axes('position', [0.32  0.2    1.1*[0.15    0.22]], FZ, 12);
        %Uidx = 214062202;
        Uidx = 214061903;
        T = JLanova2(Uidx);
        S = JLdatastruct(Uidx);
        CS = JLcycleStats(Uidx);
        Phase = linspace(0,1-1/T.NsamPerCycle,T.NsamPerCycle);
        CC = contourf(Phase, Phase, T.var2d, CS.VarSpont1*[0 0.25 0.5 1:0.5:5]);
        hold on
        PercBelow = 100*sum(T.var2d(:)<CS.VarSpont1)/numel(T.var2d)
        contour(Phase, Phase, T.var2d, CS.VarSpont1, 'color', 'w', 'linewidth', 3)
        %contour(Phase, Phase, T.var2d, CS.VarSpont1, 'color', 'w', 'linewidth', 1)
        xlabel('Contra phase (cycle)', FZ, 13);
        title('Across-beat variance', 'fontangle', 'italic');
        local_mark(ha, 'E', -0.2, 0.3, 'backgroundcolor', 'none');
        DIPratio = min(T.var2d(:))/CS.VarSpont1
    case '4F' % monaural variance
        local_fig(3);
        Xleft = 0.6;
        ha = axes('position', [Xleft  0.2    1.1*[0.15    0.22]], FZ, 12);
        set(ha, FZ,12);
        [Uidx_I Uidx_C] = deal(214062103, 214062003);
        Wi = JLwaveforms(Uidx_I);
        Wc = JLwaveforms(Uidx_C);
        SV = retrieve_dbase('D:\processed_data\JL\JLsplitvar\JLsplitvar.dbase');
        sv = SV([SV.Uidx_I]==Uidx_I);
        VARI = 100*Wi.IpsiVar/sv.SpontVarIpsi; NI = numel(VARI);
        VARC = 100*Wc.ContraVar/sv.SpontVarContra; NC = numel(VARC);
        h(1) = plot(xlim,[100 100], 'linewidth', 3, 'color', 0.7*[1 1 1]);
        xplot(xlim,[50 50], '-', 'linewidth', 1.5, 'color', 0.6*[1 1 1]);
        h(2) = xdplot(1/NI, VARI, 'linewidth', 3, ipsicolor);
        h(3) = xdplot(1/NC, VARC, 'linewidth', 3, contracolor);
        xlabel('Phase (cycle)', FZ, 13);
        ylabel('Normalized variance (%)', FZ, 13);
        hleg = legend(h, {'spont' 'ipsi' 'contra'})
        ylim([0 249]);
        100*([min(VARI) min(VARC)])
        set(ha, 'box', 'off');
        set(hleg, 'posi', [Xleft+0.06  0.35    0.0947    0.0908])
        local_mark(ha, 'F', -0.2, 0.3, 'backgroundcolor', 'none');
        % pop
        ICbothDip = [SV.minIvar]<=[SV.SpontVarIpsi] & [SV.minCvar]<=[SV.SpontVarContra];
        percentage_dipping  =100*sum(ICbothDip)/numel(SV)
        svx = SV(ICbothDip); 
        NtotCell = numel(unique([SV.icell_run]));
        NdipCell = numel(unique([svx.icell_run]));
    case '5A' % waveforms I,C,B showing thresholding
        local_fig(5);
        ha = axes('position', [0.0300    0.5436    0.2250    0.3864], FZ, 12); 
        JL_wcartoon(214061805,8, [1514 1000 1830]-2.5);
        set(gca,'visib', 'off');
        xlim([-2 7]); ylim([-1.4 3.3]);
        %xlim([2.1059    8.9809]); ylim([-1.3155    3.2923]);
        %xplot(7.0, 1.2014, 'or', 'markersize', 8)
        text(-1, 2.7284, 'I', 'fontsize', 14, 'color', [0 0 0.8], 'horizontalalign', 'center', 'fontweight', 'bold');
        text(-1, 1.7284, 'C', 'fontsize', 14, 'color', [0.8 0 0 ], 'horizontalalign', 'center', 'fontweight', 'bold');
        text(-1, 0, 'BIN', 'fontsize', 14, 'color', [0 0 0], 'horizontalalign', 'center', 'fontweight', 'bold');
        text(5.0, 1.1434, '\downarrow', 'fontsize', 17, 'fontweight', 'bold' ,'rotation', -25)
        % [hb, htx, hty] = local_scalehook(XY, Lw, DX, tX, dxt, DY, tY, dyt, FZ);
        [hb, htx, hty] = local_scalehook([2.8 -1], 3, -2, '2 ms', 0.1, 0.5, '0.5 mV', 0.2, 12);
        local_mark(ha, 'A', 0, 0.205, 'backgroundcolor', 'none');
    case '5B' % hist of ratio binaural/(summed monaural) spike rates
        local_fig(5);
        ha = axes('position', [0.3312    0.6000    0.2250    0.3300], FZ, 12); 
        C = local_triplets();
        Y = C;
        N = numel(C);
        MonRate = [C.DrivenRate_I]+[C.DrivenRate_C]-[C.SpontRate1_I]-[C.SpontRate1_C];
        BinRate = [C.Rate_B_best]-[C.SpontRate1_B];
        %BinRate = [C.Rate_B_zero]-[C.SpontRate1_B];
        Rat = BinRate./MonRate;
        Rat = Rat(Rat>0);
        local_loghist(Rat, 11, [0.5 1 2 5 10 20]);
        xlim(log([0.2 100]));
        ylim([0 24]);
        hl=findobj(ha,'type', 'line', 'marker', '*'); delete(hl);
        MM = geomean(Rat);
        set(ha, FZ, 12, 'box', 'off');
        xlabel('Binaural/monaural spike count ratio', FZ, 14);
        ylabel('# per bin', FZ, 14);
        text(0.58, 0.91, sprintf('Mean = %1.1f', MM), 'units', 'normalized', FZ, 13, ...
            'backgroundcolor', 'w');
        local_mark(ha, 'B', -0.1, 0.2, 'backgroundcolor', 'none');
    case '5C',% power-law example
        local_fig(5);
        ha = axes('position', [0.6606    0.6000    0.2250    0.3300], FZ, 12); 
        %Uidx_B = 214062808; % Exp 214 cell 6  --- B  800 Hz   40 dB (series 6)
        %Uidx_B = 214062505; % Exp 214 cell 6  --- B  500 Hz   50 dB (series 7)
        %Uidx_B = 214061906; % Exp 214 cell 6  --- B  600 Hz   70 dB  (series 9)
        %Uidx_B = 214062505; % Exp 214 cell 6  --- B  500 Hz   50 dB  (series 7)
        Uidx_B = 214062203;% Exp 214 cell 6  --- B  300 Hz   80 dB  (series 10)
        SV = retrieve_dbase('D:\processed_data\JL\JLsplitvar\JLsplitvar.dbase');
        sv = SV([SV.Uidx_B]==Uidx_B);
        Uidx_I = sv.Uidx_I; Uidx_C = sv.Uidx_C;
        P_B = JLpowerlaw(Uidx_B);
        P_I = JLpowerlaw(Uidx_I);
        P_C = JLpowerlaw(Uidx_C);
        plot(P_B.Vinst-P_B.MeanDriv, P_B.cond_spike_prob, 'o-', 'color', 'k', MFC, 'k', 'markersize', 8);
        xplot(P_C.Vinst-P_C.MeanDriv, P_C.cond_spike_prob, '^-', contracolor, 'markersize', 8);
        xplot(P_I.Vinst-P_I.MeanDriv, P_I.cond_spike_prob, 'v-', ipsicolor, 'markersize', 8);
        legend({'Binaural' 'Ipsi' 'Contra'}, 'location', 'northwest', FZ, 12);
        set(ha, FZ, 12, 'box', 'off');
        xlabel('Potential (mV)', FZ, 14);
        ylabel('Firing probability (sp/s)', FZ, 14);
        local_mark(ha, 'C', -0.1, 0.2, 'backgroundcolor', 'none');
    case 'X4D',% power-law example
        % 223/4 500 Hz, 50 dB 
        BIC = [223041204 223041404 223041304]+1;
        Wb = JLwaveforms(BIC(1)); CSb = JLcycleStats(BIC(1));
        Wi = JLwaveforms(BIC(2)); CSi = JLcycleStats(BIC(2));
        Wc = JLwaveforms(BIC(3)); CSc = JLcycleStats(BIC(3));
        figure; set(gcf,'units', 'normalized', 'position', [0.28 0.215 0.315 0.685]);
        %==B==
        ha(1)=subplot(3,1,3);
        mm = Wb.BeatMeanrec; mm = mm-mean(mm); sd = sqrt(Wb.BeatVar); dt = Wb.dt_BeatMean;
        %dplot(dt, mm); xdplot(dt, mm+sd,'r'); xdplot(dt, mm-sd,'r'); 
        local_featherplot(dt, mm, sd, [0 0.6 0], [0 0 0]);
        xlim([138.5 142.4]); 
        xlim(139.6+[0 2]); 
        DT = diff(xlim); T0 = min(xlim); T1 = max(xlim);
        %==I==
        ha(2)=subplot(3,1,1);
        mm = Wi.IpsiMeanrec; mm = mm-mean(mm); sd = sqrt(Wi.IpsiVar); dt = [Wi.dt_IpsiMean T0];
        mm = repmat(mm,[5 1]); sd = repmat(sd,[5 1]); 
        %dplot(dt, mm); xdplot(dt, mm+sd,'r'); xdplot(dt, mm-sd,'r'); 
        local_featherplot(dt, mm, sd, [0 0 0.6], [0 0 0.8]);
        DDT = 1.4;
        xlim([T0 T1]+2.73);
        %==C==
        ha(3)=subplot(3,1,2);
        mm = Wc.ContraMeanrec; mm = mm-mean(mm); sd = sqrt(Wc.ContraVar); dt = [Wc.dt_ContraMean T0];
        mm = repmat(mm,[5 1]); sd = repmat(sd,[5 1]); 
        %dplot(dt, mm); xdplot(dt, mm+sd,'r'); xdplot(dt, mm-sd,'r'); 
        local_featherplot(dt, mm, sd, [0.6 0 0], [0.8 0 0]);
        xlim([T0 T1]+1.6);
        for ii=1:3, ylim(ha(ii), [-0.4 0.65]); end
        set(ha,'xtick', [])
        text(0.4, 0.9, 'BIN', 'units', 'normalized', 'fontsize', 14, 'parent', ha(1))
        text(0.4, 0.9, 'IPSI', 'units', 'normalized', 'fontsize', 14, 'parent', ha(2))
        text(0.4, 0.9, 'CONTRA', 'units', 'normalized', 'fontsize', 14, 'parent', ha(3))
        set(gcf,'PaperPositionMode', 'auto')
    case '4C_old' % predicting monaural spike rates from binaural power laws
        local_fig(4);
        ha = axes('position', [0.7  0.6    1.5*[0.15    0.22]], FZ, 12); 
        PR = local_monspikepred();
        Y = PR;
        N = numel(PR);
        plot([PR.Nsp_I], [PR.Nsp_pred_I], '.', 'color', [0.2 0.2 1], 'markersize', 9); 
        xplot([PR.Nsp_C], [PR.Nsp_pred_C], '.', contracolor, 'markersize', 10);
        xlim([0 400]); ylim(xlim);
        xlog125([5 1200]);
        ylog125([5 1200]);
        xlabel('# spikes observed', FZ, 14);
        ylabel('# spikes predicted', FZ, 14);
        CorI = corr([PR.Nsp_I]', [PR.Nsp_pred_I]');
        CorC = corr([PR.Nsp_C]', [PR.Nsp_pred_C]');
        %text(0.2, 0.8, sprintf('R=(%0.2f,%0.2f)', CorI, CorC), FZ, 13, 'units', 'normalized');
        legend({sprintf('Ipsi R=%0.2f', CorI), sprintf('Contra R=%0.2f', CorC)}, 'location', 'northwest', FZ, 12);
        xplot(xlim, xlim, 'k')
        local_mark(ha, 'C', -0.1, 0.2, 'backgroundcolor', 'none');
    case '4X' % predicting monaural spike rates from binaural power laws
        figure;
        set(gcf,'units', 'normalized', 'position', [0.28 0.236 0.595 0.663]);
        ha = gca;
        PR = local_monspikepred();
        Y = PR;
        N = numel(PR);
        h_i = plot([PR.Nsp_I], [PR.Nsp_pred_I], '.', 'color', [0.2 0.2 1], 'markersize', 9); 
        h_c = xplot([PR.Nsp_C], [PR.Nsp_pred_C], '.', contracolor, 'markersize', 10);
        xlim([0 400]); ylim(xlim);
        %xlog125([5 1200]);
        %ylog125([5 1200]);
        xlabel('# spikes observed', FZ, 14);
        ylabel('# spikes predicted', FZ, 14);
        CorI = corr([PR.Nsp_I]', [PR.Nsp_pred_I]');
        CorC = corr([PR.Nsp_C]', [PR.Nsp_pred_C]');
        %text(0.2, 0.8, sprintf('R=(%0.2f,%0.2f)', CorI, CorC), FZ, 13, 'units', 'normalized');
        legend({sprintf('Ipsi R=%0.2f', CorI), sprintf('Contra R=%0.2f', CorC)}, 'location', 'northwest', FZ, 12);
        xplot(xlim, xlim, 'k')
        local_mark(ha, 'C', -0.1, 0.2, 'backgroundcolor', 'none');
        IDpoints(h_i, 'IPSI', PR, @(s,pr)sprintf('%s  %s', s, JLtitle(pr.Uidx_B)), 'powerfun', @local_powerfun);
    case '4D', % histograms of pred/observed Nspikes
        set(figure,'units', 'normalized', 'position', [0.28 0.226 0.63 0.674]);
        PR = local_monspikepred();
        Y = PR;
        rat_I = deinf(A2dB([PR.Nsp_pred_I]./[PR.Nsp_I]));
        rat_I = rat_I(abs(rat_I)<399);
        rat_C = deinf(A2dB([PR.Nsp_pred_C]./[PR.Nsp_C]));
        rat_C = rat_C(abs(rat_C)<399);
        subplot(3,2,1); 
        hist(rat_I,20); 
        xlabel('Ipsi #spikes pred/obs (dB)')
        text(0.2, 0.9, sprintf('median ratio %0.2f', dB2A(nanmedian(rat_I))), 'units', 'normali');
        subplot(3,2,2); 
        hist(rat_C,20);
        set(findobj(gca, 'type', 'patch'), 'FaceColor', 'r');
        xlabel('Contra #spikes pred/obs (dB)')
        min(rat_C)
        text(0.2, 0.9, sprintf('median ratio %0.2f', dB2A(nanmedian(rat_C))), 'units', 'normali');
        subplot(3,2,3); 
        plot([PR.Nsp_I], [PR.Nsp_C], 'k.');
        xlog125([10 1000]); ylog125([10 1000]);
        xlabel('Obs #spikes ipsi');
        ylabel('Obs #spikes contra');
        %axis equal
        xplot(xlim, xlim, 'k')
        subplot(3,2,4); 
        plot([PR.Nsp_pred_I], [PR.Nsp_pred_C], 'k.');
        xlog125([10 1000]); ylog125([10 1000]);
        xlabel('Pred #spikes ipsi');
        ylabel('Pred #spikes contra');
        %axis equal
        xplot(xlim, xlim, 'k')
        %subplot(3,2,5); 
    case '102' % best ITD varies w freq
        local_fig(102);
        PS = struct('linestyle', '-', 'marker', 'o', 'color', [0 0 0.7], 'markerfacecolor', [0 0 0.7], 'linewidth',2, 'markersize', 8)
        AT = @(B)text(0.95, 0.05, local_RG(B), 'units', 'normalized', 'fontsize', 13, 'color', 0.7*[1 1 1], 'horizontalalign', 'right');
        LM = @(ha, CH)local_mark(ha, CH, -0.1, 0.2, 'backgroundcolor', 'none');
        %
        ha(1) = subplot(2,2,1);
        B = JLbestITD(214061800+(1:15));
        plot([B.freq], [B.bestITD], PS);
        xlim([0 1000]);
        ylim([-500 1000]);
        ylabel('Best ITD (\mus)', FZ,15);
        AT(B);
        LM(gca,'A');
        set(gca,'xtick',0:250:1250);
        xplot(xlim', 130*[-1 1; -1 1], '--', 'color', [0 0.8 0], 'linewidth', 3);
        %
        ha(2) = subplot(2,2,2);
        B = JLbestITD(204113100+(1:12));
        plot([B.freq], [B.bestITD], PS);
        xlim([0 1300]);
        ylim([-1000 500]);
        AT(B);
        LM(gca,'B');
        set(gca,'xtick',0:250:1250);
        xplot(xlim', 130*[-1 1; -1 1], '--', 'color', [0 0.8 0], 'linewidth', 3);
        %
        ha(3) = subplot(2,2,3);
        B = JLbestITD(209051700+(1:7));
        plot([B.freq], [B.bestITD], PS);
        xlim([0 1000]);
        ylim([-500 500]);
        xlabel('Ipsi frequency (Hz)', FZ,15);
        ylabel('Best ITD (\mus)', FZ,15);
        AT(B);
        LM(gca,'C');
        set(gca,'xtick',0:250:1250);
        xplot(xlim', 130*[-1 1; -1 1], '--', 'color', [0 0.8 0], 'linewidth', 3);
        %
        ha(4) = subplot(2,2,4);
        B = JLbestITD(219065300+(1:12))
        plot([B.freq], [B.bestITD], PS);
        xlim([0 1000]);
        ylim([0 1000]);
        xlabel('Ipsi Frequency (Hz)', FZ,15);
        AT(B);
        LM(gca,'D');
        set(gca,'xtick',0:250:1250);
        xplot(xlim', 130*[-1 1; -1 1], '--', 'color', [0 0.8 0], 'linewidth', 3);
        %
        set(ha,'fontsize',13);
        saveas(gcf, '\\Ee1285a\MSO$\Figure_102.fig');
    case '103' % ITD population data
        local_fig(103);
        ha = axes('position', [0.65  0.3    1.1*[0.15    0.22]]);
        qqq = local_IPDstuff;
        qqq = qqq(ismember([qqq.UniqueRecordingIndex], Usbu));
        % reduce qqq to best-spiking members of each cell
        icr = JL_icell_run([qqq.UniqueRecordingIndex]);
        uicr = unique(icr);
        %         QQQ = [];
        %         for ii=1:numel(uicr),
        %             qq_ii = qqq(icr==uicr(ii));
        %             [dum, isort] = sort([qq_ii.Nspikes]);
        %             i0 = max(1,numel(isort)-0);
        %             QQQ = [QQQ qq_ii(isort(i0:end))];
        %         end
        %         numel(QQQ);
        %         qqq = QQQ;
        bestITD = 1e6*[qqq.measBestIPD]./[qqq.Freq];
        hist(bestITD,linspace(-1100, 1100, 31));
        xlim([-1100 1100]);
        fenceplot(130*[-1 1],ylim, '--', 'color', [0 0.8 0], 'linewidth', 3);
        set(ha,  FZ, 13, 'box', 'off', 'xtick', [-1000:1000:1000]);
        ylabel('Count per bin', FZ, 15);
        xlabel('best ITD (\mus)', FZ, 15);
        %local_mark(ha, 'F', -0.2, 0.3, 'backgroundcolor', 'none');
        PercentOutside = 100*sum(abs(bestITD)>130)/numel(bestITD)
        MED_ITD = median(bestITD)
        IQR_ITD = iqr(bestITD)
        MEAN_ITD = mean(bestITD)
        STD_ITD = std(bestITD)
        saveas(gcf, '\\Ee1285a\MSO$\Figure_103.fig');
    case '104', % monaural & binaural phase/freq contours
        local_fig(104);
        ha = axes('position', [0.42  0.27    0.22    0.38]);
        [S, mWvI, mWvC] = JLpeakPhase(-214062507, 5.5);
        local_PPH_plot(S, mWvI, mWvC);
        saveas(gcf, '\\Ee1285a\MSO$\Figure_104.fig');
    case '105', % SBC cycle histograms
        uiopen('\\Ee1285a\MSO$\Thomas_cyclehist_ap.fig',1);
        figh = gcf;
        set(figh,'tag', 'MSO_fig105', 'units', 'pixels', 'position', [94  213  1066 734], ...
            'color', 'w', 'defaultaxestickdir', 'out', ...
            'PaperPositionMode', 'auto', 'papersize', [11 8.5]);
        hax = findobj(figh, 'type', 'axes');
        hax = flipud(hax);
        delete(hax(11:end));
        freq = [156 186 221 263 313 372 442 526 625 743];
        for ii=6:10,
            Pos = get(hax(ii), 'posi');
            Pos([2 3 4]) = [0.2 0.1 0.2];
            set(hax(ii) ,'position', Pos, 'linewidth', 2, 'box', 'off', FZ, 13);
        end
        for ii=1:5,
            Pos = get(hax(ii), 'posi');
            Pos([2 3 4]) = [0.6 0.1 0.2];
            set(hax(ii) ,'position', Pos, 'linewidth', 2, 'box', 'off', FZ, 13);
        end
        MCH = 'ABCDEFGHIJ';
        for ii=1:10,
            axes(hax(ii));
            title('');
            text(0.1, 18, sprintf('%d Hz', freq(ii)), 'fontsize', 13);
            local_mark(gca, MCH(ii), 0.1, 0.2, 'backgroundcolor', 'none');
            XL = ''; YL = '';
            if ii>5, XL = 'Phase (cycle)'; end
            if rem(ii,5)==1, YL = '#APs'; end
            xlabel(XL, FZ, 15)
            ylabel(YL, FZ, 15)
        end
        saveas(gcf, '\\Ee1285a\MSO$\Figure_105.fig');
    case '106', % 3 more bin receptive fields 
        local_fig(106);
        LM = @(ha,CH)local_mark(ha, CH,0,0.12);
        subplot(1,3,1);
        Uidx = 204112801; JLtitle(Uidx)
        [S, mWvI, mWvC] = JLpeakPhase(-Uidx, 5.5);
        local_BRF_plot(S, mWvI, mWvC);
        LM(gca, 'A')
        %=
        subplot(1,3,2);
        Uidx = 214061804; JLtitle(Uidx)
        [S, mWvI, mWvC] = JLpeakPhase(-Uidx, 5.5);
        local_BRF_plot(S, mWvI, mWvC,0);
        LM(gca, 'B')
        %=
        subplot(1,3,3);
        Uidx = 219065301; JLtitle(Uidx)
        [S, mWvI, mWvC] = JLpeakPhase(-Uidx, 5.5);
        local_BRF_plot(S, mWvI, mWvC,0);
        LM(gca, 'C')
        %=
        saveas(gcf, '\\Ee1285a\MSO$\Figure_106.fig');
    case '107', % dev btw mon sum pred and bin subthr
        LM = @(ha,CH)local_mark(ha, CH,0,0.18);
        local_fig(107);
        set(gcf,'units', 'normalized', 'position', [0.0727 0.596 0.9 0.331])
        C = local_triplets();
        N = numel(C);
        Y = C;
        %saveas(gcf, '\\Ee1285a\MSO$\Figure_107.fig');
        Uidx_B = 214061903;
        Uidx_I = 214062103
        Uidx_C = 214062003
        T = JLanova2(Uidx_B);
        nshift = [0 0];%*nshift;
        BIN = circshift(T.mean2d, nshift);
        BIN = BIN-mean(BIN(:));
        T = JLanova2(Uidx_I);
        PHi = linspace(0,1-1/T.NsamPerCycle,T.NsamPerCycle);
        Wmi = circshift(T.mean2d, nshift);
        Wmi = mean(Wmi,2);
        Wmi = Wmi-mean(Wmi(:));
        T = JLanova2(Uidx_C);
        PHc = linspace(0,1-1/T.NsamPerCycle,T.NsamPerCycle);
        nshift = [0 0];%*nshift;
        Wmc = circshift(T.mean2d, nshift);
        Wmc = mean(Wmc,1);
        Wmc = Wmc-mean(Wmc(:));
        Wmc = Wmc-mean(Wmc(:));
        Ni = numel(Wmi); PHi = (0:Ni-1)/Ni;
        Nc = numel(Wmc); PHc = (0:Nc-1)/Nc;
        [Wmc, Wmi] = SameSize(Wmc(:).', Wmi(:));
        linpred = Wmi+Wmc;
        %==
        AXpos = [
            0.1395    0.1672    0.2020    0.73
            0.4203    0.1672    0.2020    0.73
            0.7011    0.1672    0.2020    0.73];
        %=
        axes('position', AXpos(1,:));
        XL = [-0.5 1]+3;
        plot(XL,XL,'k','linewidth',2);
        xplot(3+linpred(:), 3+BIN(:), '.', 'markersize', 10);
        set(gca,FZ,13);
        xlabel('Linear prediction (mV)', FZ, 15)
        ylabel('Binaural recording (mV)', FZ, 15)
        text(0.1,0.9, sprintf('Var_{acc}=%2.1f %%', 95.5), 'units', 'normal', FZ, 13);
        LM(gca,'A');
        %=
        axes('position', AXpos(2,:));
        DEV = BIN(:)-linpred(:); 
        hist(DEV,20);
        SK = skewness(DEV);
        xlim([-0.2 0.2]);
        set(gca,FZ,13);
        xlabel('Prediction - recording (mV)', FZ, 15)
        ylabel('Count per bin', FZ, 15)
        LM(gca,'B');
        %=
        axes('position', AXpos(3,:));
        contourf(PHc, PHi, BIN-linpred ,-0.2:0.05:0.2);
        set(gca,FZ,13);
        xlabel('Contra Phase (cycle)');
        ylabel('Ipsi Phase (cycle)');
        LM(gca,'C');
        saveas(gcf, '\\Ee1285a\MSO$\Figure_107.fig');
    case '108', % test of Jercog's asymmetry hypothesis
        local_fig(108);
        ha = subplot(2,2,1);
        DBFN = 'D:\processed_data\JL\Jercog\JL_jercog.dbase';
        qq = retrieve_dbase(DBFN);
        qqc = qq([qq.P2PMeanIpsiCycle]>2*sqrt([qq.VarSpont1]) & [qq.P2PMeanContraCycle]>2*sqrt([qq.VarSpont1]));
        AllCells = unique([qqc.icell_run]);
        plot([0 4], [0 4], 'color', 0.8*[1 1 1], 'linewidth',2)
        for ii=1:numel(AllCells),
            qqc_i = qqc([qqc.icell_run]==AllCells(ii));
            xplot([qqc_i.maxSlope_I], [qqc_i.maxSlope_C], [ploco(ii) ploma(ii)], 'linewidth',2, 'markersize', 6);
        end
        set(gca,FZ, 13);
        xlim([0 3.7]); ylim([xlim]);
        axis equal; axis square;
        xlabel('Ipsi slope (V/s)', FZ, 13);
        ylabel('Contra slope (V/s)', FZ, 13);
        MeanDiff = mean([qqc.maxSlope_I]-[qqc.maxSlope_C])
        std([qqc.maxSlope_I]-[qqc.maxSlope_C])
        [TT, p] = ttest([qqc.maxSlope_I]-[qqc.maxSlope_C])
        set(gca,'box', 'off')
        MDY = 0.1;
        local_mark(ha, 'A', 0, MDY);
%    case '109' % testing Zhou's asymmetry
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
        %         fh = figure;
        %         set(fh,'units', 'pixels', 'position', [255 348 388 349], ...
        %             'color', 'w', 'defaultaxestickdir', 'out', ...
        %             'PaperPositionMode', 'auto', 'papersize', [11 8.5]);
        ha = subplot(2,2,2);
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
            xplot(meanlat_I(ii), meanlat_C(ii), [ploco(ii) ploma(ii)], 'markersize', 8, 'linewidth', 2);
        end
        TT = ttest(denan(meanlat_I-meanlat_C))
        axis equal; axis square;
        xlim([0 388]); ylim(xlim);
        xplot(xlim,xlim, 'k');
        Y = CollectInStruct(meanlat_I, meanlat_C);
        set(gca, 'fontsize', 12)
        xlabel('Ipsi EPSP-AP latency (\mus)', 'fontsize' ,13)
        ylabel('Contra EPSP-AP latency (\mus)', 'fontsize' ,13)
        local_mark(ha, 'B', 0, MDY);
    case '110' % below-spont dip in beat-cycle-avg var. More of 3E; pop plot (e.g. hist of dip/spont or dip/spont vs freq).
        local_fig(110);
        214041005; 
        %Uidx = [214061806; 216020804; 219010107; 219031902]; 
        %Uidx = [219055102; 219065303; 223020407; 223030803]; 
        %Uidx = [223041205; 223041205; 223041205; 223041205];
        Uidx = [214062202 219055102 219065303 204125802 204136703 209051703 214041005 ];
        PhSh = {[0 0] [0 0] [0 0.5] [0.5 0.5]};
        %==
        DX = 0.22;
        [DIPratio(1), ha] = local_binvardip(Uidx(1), 'AE', 0.1, PhSh{1});
        xlabel(ha(2), 'Contra phase (cycle)', FZ, 13);
        ylabel(ha(1), 'Ipsi phase (cycle)', FZ, 13);
        ylabel(ha(2), 'Ipsi phase (cycle)', FZ, 13);
        %==
        [DIPratio(2), ha] = local_binvardip(Uidx(2), 'BF', 0.1+DX, PhSh{2});
        xlabel(ha(2), 'Contra phase (cycle)', FZ, 13);
        %==
        [DIPratio(3), ha] = local_binvardip(Uidx(3), 'CG', 0.1+2*DX, PhSh{3});
        xlabel(ha(2), 'Contra phase (cycle)', FZ, 13);
        %==
        [DIPratio(4), ha] = local_binvardip(Uidx(4), 'DH', 0.1+3*DX, PhSh{4});
        xlabel(ha(2), 'Contra phase (cycle)', FZ, 13);
        DIPratio
    case '111' % time spent below spont var as a fnc of freq
        local_fig(111);
        UidxBase = 214062200; % 80 dB
        %UidxBase = 214062500; % 50 dB
        %UidxBase = 214061900; % 70 dB
        for ii=1:15,
            Uidx=UidxBase+ii;
            CS = JLcycleStats(Uidx);
            SPv(ii) = CS.VarSpont1;
        end
        SPv = mean(SPv);
        for ii=1:15,
            Uidx=UidxBase+ii;
            T = JLanova2(Uidx);
            S = JLdatastruct(Uidx);
            Freq(ii) = S.Freq1;
            PercBelow(ii) = 100*sum(T.var2d(:)<SPv)/numel(T.var2d);
        end
        ha = subplot(2,2,1);
        set(ha, FZ, 14);
        xplot(Freq, PercBelow, 'k*-', 'linewidth', 2)
        xlabel('Frequency (Hz)', FZ, 17);
        %ylabel(['Fraction of beat cycle' char(10) 'having \sigma^2<\sigma^2_{spont} (%)'], FZ, 13);
        ylabel('Fraction of beat cycle (%)', FZ, 17);
        local_mark(ha, 'A', -0.1, 0.2, 'backgroundcolor', 'none', 'fontsize', 20);
        UidxBase = 214062400; % 80 dB, I
        for ii=1:15,
            Uidx=UidxBase+ii;
            CS = JLcycleStats(Uidx);
            SPv(ii) = CS.VarSpont1;
        end
        SPv = mean(SPv);
        for ii=1:15,
            Uidx=UidxBase+ii;
            T = JLanova2(Uidx);
            S = JLdatastruct(Uidx);
            Freq(ii) = S.Freq1;
            PercBelow(ii) = 100*sum(T.var2d(:)<SPv)/numel(T.var2d);
        end
        ha = subplot(2,2,2);
        set(ha, FZ, 14);
        xplot(Freq, PercBelow, 'k*-', 'linewidth', 2)
        xlabel('Frequency (Hz)', FZ, 17);
        %ylabel(['Fraction of beat cycle' char(10) 'having \sigma^2<\sigma^2_{spont} (%)'], FZ, 13);
        ylabel('Fraction of beat cycle (%)', FZ, 17);
        local_mark(ha, 'B', -0.1, 0.2, 'backgroundcolor', 'none', 'fontsize', 20);
        %UidxBase = 214062300; % 80 dB, C
    case 'S3B', % size of cycle-averaged waveforms. Mon versus Bin 
        figure;
        ha = subplot(2,2,1);
        SV = retrieve_dbase('D:\processed_data\JL\JLsplitvar\JLsplitvar.dbase');
        SV = SV(abs(A2dB([SV.MeanSpontVar]./[SV.SpontVarBin]))<4);
        plot([SV.SPL]+1.5*randn(size(SV)), ([SV.RMSbinIpsi]./sqrt([SV.SpontVarBin]))./([SV.RMSmonIpsi]./sqrt([SV.SpontVarIpsi])), '.')
        xplot([SV.SPL]+1.5*randn(size(SV)), ([SV.RMSbinContra]./sqrt([SV.SpontVarBin]))./([SV.RMSmonContra]./sqrt([SV.SpontVarContra])), 'r.')
        ylog125([0.1 10]);
        xplot(xlim, [1 1], 'k');
        set(ha, FZ, 12)
        xlabel('Intensity (dB SPL)', FZ, 13);
        ylabel('Bin/Mon amp ratio', FZ, 13);
        local_mark(ha, 'S3B', -0.1, 0.15, 'backgroundcolor', 'none');
    case 'S3C' % pop data for reduction of binaural var below spont
        SV = retrieve_dbase('D:\processed_data\JL\JLsplitvar\JLsplitvar.dbase');
        SV = SV([SV.freq]<1800);
        SV = SV(rem([SV.freq],100)==0);
        SV30 = SV([SV.SPL]==30); SV40 = SV([SV.SPL]==40); SV50 = SV([SV.SPL]==50);
        SV60 = SV([SV.SPL]==60); SV70 = SV([SV.SPL]==70); SV80 = SV([SV.SPL]==80);
        cla;
        ha = subplot(2,2,1);
        %xplot([SV30.freq], 100*[SV30.minBinVar]./[SV30.SpontVarBin],'*', lico(1));
        %xplot([SV40.freq], 100*[SV40.minBinVar]./[SV40.SpontVarBin],'*', lico(2));
        %xplot([SV50.freq], 100*[SV50.minBinVar]./[SV50.SpontVarBin],'*', lico(3));
        xplot([SV60.freq], 100*[SV60.minBinVar]./[SV60.SpontVarBin],'*', lico(4));
        xplot([SV70.freq], 100*[SV70.minBinVar]./[SV70.SpontVarBin],'*', lico(5));
        xplot([SV80.freq], 100*[SV80.minBinVar]./[SV80.SpontVarBin],'*', lico(6));
        xlim([0 1300])
        set(ha, FZ, 12)
        xlabel('Frequency (Hz)', FZ, 13);
        ylabel('BinDip var over Spont var (%)', FZ, 13);
        hleg = legend({'30 dB', '40 dB', '50 dB', '60 dB', '70 dB', '80 dB'}, 'location', 'northeast')
        set(hleg, FZ, 10)
        local_mark(ha, 'S3C', -0.1, 0.15, 'backgroundcolor', 'none');
    case 'S3D' % cell 30: reduction of dip monaural var re spont var against Cycle-mean size
        SV = retrieve_dbase('D:\processed_data\JL\JLsplitvar\JLsplitvar.dbase');
        sv = SV([SV.icell_run]==30);
        set(figure,'units', 'normalized', 'position', [0.28 0.256 0.26 0.644])
        ha1 = subplot(2,1,1); 
        plot([sv.MeanSpontVar]-[sv.minIvar], [sv.RMSmonIpsi], '.'); 
        xplot([sv.MeanSpontVar]-[sv.minCvar], [sv.RMSmonContra], 'r.'); 
        xlabel('Mon Dip var reduction re spont (mV^2)');
        legend({'Ipsi' 'Contra'});
        ylabel('RMS mon cycle (mV)');
        local_mark(ha1, 'S3D1', -0.1, 0.15, 'backgroundcolor', 'none');
        ha2 = subplot(2,1,2); 
        plot([sv.MeanSpontVar]-[sv.minIvar], [sv.MeanSpontVar]-[sv.minCvar], '.');
        xlabel('Ipsi Dip var reduction re spont (mV^2)');
        ylabel('Contra Dip var reduction re spont (mV^2)');
        xplot(xlim, xlim, 'k'); corr([sv.MeanSpontVar]'-[sv.minIvar]', [sv.RMSmonIpsi]')
        axis square;
        local_mark(ha2, 'S3D2', -0.1, 0.15, 'backgroundcolor', 'none');
        set(gcf, 'PaperPositionMode', 'auto')
    otherwise
        error('unknown figure panel.');
end


%====================
function h = local_fig(k);
FN = ['MSO_fig' num2str(k)];;
switch k
    case 1,
        h = findobj('type', 'figure', 'tag', FN);
        if isempty(h),
            h = figure;
            set(h,'tag', FN, 'units', 'pixels', 'position', [94  213  1066 734], ...
                'color', 'w', 'defaultaxestickdir', 'out', ...
                'PaperPositionMode', 'auto', 'papersize', [11 8.5]);
        end
    case 2,
        h = findobj('type', 'figure', 'tag', FN);
        if isempty(h),
            h = figure;
            set(h,'tag', FN, 'units', 'pixels', 'position', [94  213  1066 734], ...
                'color', 'w', 'defaultaxestickdir', 'out', ...
                'PaperPositionMode', 'auto', 'papersize', [11 8.5]);
        end
    case 3,
        h = findobj('type', 'figure', 'tag', FN);
        if isempty(h),
            h = figure;
            set(h,'tag', FN, 'units', 'pixels', 'position', [94  213  1066 734], ...
                'color', 'w', 'defaultaxestickdir', 'out', ...
                'PaperPositionMode', 'auto', 'papersize', [11 8.5]);
        end
    case 4,
        h = findobj('type', 'figure', 'tag', FN);
        if isempty(h),
            h = figure;
            set(h,'tag', FN, 'units', 'pixels', 'position', [94  213  1066 734], ...
                'color', 'w', 'defaultaxestickdir', 'out', ...
                'PaperPositionMode', 'auto', 'papersize', [11 8.5]);
        end
    case 5,
        h = findobj('type', 'figure', 'tag', FN);
        if isempty(h),
            h = figure;
            set(h,'tag', FN, 'units', 'pixels', 'position', [94  213  1066 734], ...
                'color', 'w', 'defaultaxestickdir', 'out', ...
                'PaperPositionMode', 'auto', 'papersize', [11 8.5]);
        end
    case {102,103,104,105,106,107,108,109,110,111,112},
        h = findobj('type', 'figure', 'tag', FN);
        if isempty(h),
            h = figure;
            set(h,'tag', FN, 'units', 'pixels', 'position', [94  213  1066 734], ...
                'color', 'w', 'defaultaxestickdir', 'out', ...
                'PaperPositionMode', 'auto', 'papersize', [11 8.5]);
        end
    otherwise, error('invalid fig number')
end
figure(h);

function P = local_panmarker(ha);
P = struct('parent', ha, 'units', 'normalized', ...
    'horizontalalign', 'center', ...
    'fontsize', 13 ,'fontweight', 'bold' ,'BackgroundColor', 'w');

function local_mark(ha, C, DX, DY, varargin);
if nargin<3, DX=0; end
if nargin<4, DY=0; end
text(0.04+DX, 0.9+DY, lower(C), local_panmarker(ha), varargin{:});

function [hb, ht] = local_scalebar(X,Y, DX, DY, Lw ,Txt, dyt, Fsz);
hb = xplot(X+DX/2*[-1 -1 -1 1 1 1], Y+DY/2*[1 -1 0 0 1 -1], 'k', 'linewidth', Lw);
ht = text(X, Y+dyt, Txt, 'fontsize', Fsz, 'fontweight', 'bold', 'horizontalalign', 'center');

function [hb, ht] = local_scalebarY(X,Y, DX, DY, Lw ,Txt, dxt, Fsz);
hb = xplot(X+DX/2*[1 -1 0 0 1 -1], Y+DY/2*[-1 -1 -1 1 1 1], 'k', 'linewidth', Lw);
ht = text(X+dxt, Y, Txt, 'fontsize', Fsz, 'fontweight', 'bold', 'horizontalalign', 'left');

function [hb, htx, hty] = local_scalehook(XY, Lw, DX, tX, dxt, DY, tY, dyt, FZ);
%XY, Lw, DX, tX, dxt, DY, tY, dyt, FZ
[X,Y] = DealElements(XY);
hb = xplot(X+[DX 0 0], Y+[0 0 DY], 'k', 'linewidth', Lw);
htx = text(X+DX/2, Y-sign(DY)*dxt, tX, 'horizontalalign', 'center', 'verticalalign' ,'top', 'fontsize', FZ);
if DX>0, HorAl = 'right'; else, HorAl = 'left'; end;
hty = text(X-sign(DX)*dyt, Y+DY/2, tY, 'horizontalalign', HorAl, 'fontsize', FZ);

function local_APthrplot(Y);
Uidx = Y.UniqueRecordingIndex;
% ==histogram of slopes, working exclusively from raw slope stats in variable minSlp
logSlope = log(Y.minSlp(Y.minSlp>0));
Nbin = min(50, numel(logSlope)/10);
BC = linspace(min(logSlope), max(logSlope), Nbin);
hist(logSlope, BC); NN = hist(logSlope, BC); % 1st call: force plot; 2nd: get hist values.
xlim([min(BC) max([max(BC), 1.05*log(Y.SlopeThr)])]);
Xval = [1 2 3 4 5 7 8 10 12 15 20 25 30 50 100];
set(gca,'xtick', log(Xval), 'xticklabel', Words2cell(num2str(Xval)), 'fontsize',11);
fenceplot(log(Y.SlopeThr), ylim, 'color', [0 0.75 0], 'linewidth',2);
Peak2 = max(NN(BC>=log(Y.SlopeThr)));
Ymax = min(max(NN), 10*Peak2);
ylim([0 Ymax]);
hold on;
hist(logSlope(logSlope>log(Y.SlopeThr)), BC);
hp = findobj(gca, 'type', 'patch');
set(hp(1), 'FaceColor', 'r')
%get(hp);
xlabel('Downward slope (V/s)', 'fontsize', 14);
ylabel('# occurrences', 'fontsize', 14);

%=====================================================
function [mWi, Wsel, dti] = local_cyclesnips(W, M, T, Nsel, DY, tshift)
% return cycle-avg waveform & a random selection of waterfall-shifted
% cycles of the recording
dti = T/M; % new sample period
Tmax = W.dt*numel(W.stimrec_detrended);
Ncycle = floor(Tmax/T);
Wi = interp1(Xaxis(W.stimrec_detrended, W.dt),W.stimrec_detrended, tshift+(0:M*Ncycle-1)*dti);
Wi = reshape(Wi, [M,numel(Wi)/M]);
mWi = nanmean(Wi,2);
mm = mean(mWi);
mWi = mWi - mm;
Wi = Wi-mm;
Wi = Wi(:,~any(isnan(Wi),1));
Dev = std(bsxfun(@minus, Wi, mWi));
Nc = size(Wi,2);
Wsel = Wi(:, 1:Nsel);
Wsel = Wsel + SameSize(DY*(1:size(Wsel,2)), Wsel);


function local_stackplot(dti, mWi, selWi, CL, freqlabel);
[M, Nsel] = size(selWi); Tipsi = M*dti;
phi=linspace(0,2*pi,M); Stim = sin(phi);
dplot(dti, 1.1*Nsel+Stim, CL, 'linewidth', 3);
xdplot(dti, 0.7*Nsel+0*mWi, 'k:', 'linewidth', 3);
xdplot(dti, selWi, 'k', 'linewidth',2);
xdplot(dti, -1+0*mWi, 'k:', 'linewidth', 3);
xdplot(dti, -3+4*mWi, CL, 'linewidth', 3);
xlim([0 Tipsi]);
ylim([-5 1.2*Nsel]);
text(mean(xlim), 0.92*Nsel, freqlabel, CL, 'horizontalalign', 'center')
set(gca,'visib', 'off')

function S = ipsicolor;
S = struct('color', [0 0 0.8], 'markerfacecolor', [0 0 0.8]);

function S = contracolor;
S = struct('color', [0.8 0 0], 'markerfacecolor', [0.8 0 0]);

function qqq = local_IPDstuff();
DBFN = fullfile('D:\processed_data\JL\monpeax2bestITD', 'JL_monpeax2bestITD.dbase');
qqq = retrieve_dbase(DBFN);
qqq = qqq([qqq.alphaBestIPD]<=0.001);
qqq = qqq([qqq.Nspikes]>5);
qqq = qqq(~isnan([qqq.predBestIPD]) & ~isnan([qqq.measBestIPD]));

function local_cplot(DX, DY, W, varargin);
N = numel(W);
phi = (0:N-1)/N;
W = W-mean(W(:));
xplot(DX + phi, W+DY, varargin{:});

function local_loghist(Y, Nbin, Xlab);
L = log(Y);
hist(L,Nbin);
set(gca, 'xtick', log(Xlab), 'xticklabel', Words2cell(num2str(Xlab)));

function id = local_isdriven(Nspont, SpontDur, Ndriv, DrivDur, Nsigma);
Rate_spont = Nspont./SpontDur;
Rate_driv = Ndriv./DrivDur;
sigma_spont = sqrt(Nspont)./SpontDur;
sigma_driv = sqrt(Ndriv)./DrivDur;
sigma_tot = sqrt(sigma_spont.^2+sigma_driv.^2);
id = (Rate_driv-Rate_spont)>Nsigma*sigma_tot;

function PR = local_monspikepred();
% driven triplets
C = retrieve_dbase('D:\processed_data\JL\JLcompareSpikeCounts\JLcompareSpikeCounts_10.dbase');
qbin = [C.vs_alpha_B]<0.005;
C = C(qbin);
DurSpont1 = 500; % ms
Nsigma = 1.3;
qdriv_I = local_isdriven([C.NspikeInSpont1_I], DurSpont1, [C.NspikeInStim_I], [C.AnaDur], Nsigma);
qdriv_C = local_isdriven([C.NspikeInSpont1_C], DurSpont1, [C.NspikeInStim_C], [C.AnaDur], Nsigma);
qdriv_B = local_isdriven([C.NspikeInSpont1_B], DurSpont1, [C.NspikeInStim_B], [C.AnaDur], Nsigma);
C = C((qdriv_I & qdriv_C) & qdriv_B);
N = numel(C);
Uidx_B = [C.Uidx_B];
DBFN = 'D:\processed_data\JL\PredMonSpikes\JLpredictMonRate_10.dbase';
PR = retrieve_dbase(DBFN);
PR = PR(~isnan([PR.Nsp_pred_C]));
minNsp = 7;
PR = PR([PR.Nsp_B]>=minNsp & [PR.Nsp_I]>=minNsp & [PR.Nsp_C]>=minNsp);
PR = PR(ismember([PR.Uidx_B],Uidx_B));
%PR = PR([PR.VSalpha_B]<=0.001 & [PR.VSalpha_Imon]<=0.001 & [PR.VSalpha_Cmon]<=0.001);

function C = local_triplets();
C = retrieve_dbase('D:\processed_data\JL\JLcompareSpikeCounts\JLcompareSpikeCounts_10.dbase');
qbin = [C.vs_alpha_B]<0.005;
C = C(qbin);
DurSpont1 = 500; % ms
Nsigma = 1;
qdriv_I = local_isdriven([C.NspikeInSpont1_I], DurSpont1, [C.NspikeInStim_I], [C.AnaDur], Nsigma);
qdriv_C = local_isdriven([C.NspikeInSpont1_C], DurSpont1, [C.NspikeInStim_C], [C.AnaDur], Nsigma);
qdriv_B = local_isdriven([C.NspikeInSpont1_B], DurSpont1, [C.NspikeInStim_B], [C.AnaDur], Nsigma);
C = C((qdriv_I | qdriv_C) & qdriv_B);
N = numel(C);



function XPN = local_RG(B);
D = JLdatastruct(B(1).UniqueRecordingIndex);
XPN = sprintf('%d/%d/%d', D.iexp, D.icell, D.iseq);

function figh = local_PPH_plot(S, mWvI, mWvC);
% phase-freq contours from JLpeakPhase
figh = gcf;
XY = [850 600];
CLR = get(figh,'color');
UIp = struct('fontsize', 12);
LM = @(ha,CH)local_mark(ha, CH,0,0.12);
subplot(1,3,1);
Nph = size(mWvI,1); % # phase values
Nh = round(Nph/2);
dph = 1/Nph;
PH = dph*(0:Nph-1).'-0.5;
[Wmin, Wmax] = minmax([mWvI(:); mWvC(:)]);
contourf(PH, [S.Freq_ipsi], circshift(mWvI ,Nh).');
caxis([Wmin, Wmax]);
ylabel('Frequency (Hz)', 'fontsize',14);
xlabel('Ipsi phase (cycles)', 'fontsize',14);
set(gca, 'xtick',-0.5:0.25:0.5, 'fontsize', 13);
xlim([-0.5 0.5])
LM(gca, 'A');
%=
subplot(1,3,3);
CC=contourf(PH, [S.Freq_ipsi], circshift(mWvC ,Nh).');
caxis([Wmin, Wmax]);
xlabel('Contra phase (cycles)', 'fontsize',14);
set(gca, 'xtick',-0.5:0.25:0.5, 'fontsize', 13, 'yticklabel', {});
xlim([-0.5 0.5])
LM(gca, 'C');
%=
subplot(1,3,2);
SAP = [S.ShiftAddPeak];
Freq = [S.Freq_ipsi];
contourf(PH, Freq, fliplr(circshift(SAP ,Nh).'));
xlabel('IPD (cycles)', 'fontsize',14);
caxis(max(SAP(:))+[2*(Wmin-Wmax) 0]);
xlim([-0.5 0.5])
set(gca, 'xtick',-0.5:0.25:0.5, 'fontsize', 13, 'yticklabel', {});
LM(gca, 'B');

function local_BRF_plot(S, mWvI, mWvC, makeYlabel);
% bin receptive field from JLpeakPhase
makeYlabel = arginDefaults('makeYlabel',1);
figh = gcf;
XY = [850 600];
CLR = get(figh,'color');
UIp = struct('fontsize', 12);
%subplot(1,3,1);
Nph = size(mWvI,1); % # phase values
Nh = round(Nph/2);
dph = 1/Nph;
PH = dph*(0:Nph-1).'-0.5;
[Wmin, Wmax] = minmax([mWvI(:); mWvC(:)]);
%contourf(PH, [S.Freq_ipsi], circshift(mWvI ,Nh).');
%caxis([Wmin, Wmax]);
%xlabel('Ipsi phase (cycles)', 'fontsize',14);
%set(gca, 'xtick',-0.5:0.25:0.5, 'fontsize', 13);
%xlim([-0.5 0.5])
%LM(gca, 'A');
%=
%subplot(1,3,3);
%CC=contourf(PH, [S.Freq_ipsi], circshift(mWvC ,Nh).');
%caxis([Wmin, Wmax]);
%xlabel('Contra phase (cycles)', 'fontsize',14);
%set(gca, 'xtick',-0.5:0.25:0.5, 'fontsize', 13, 'yticklabel', {});
%xlim([-0.5 0.5])
%LM(gca, 'C');
%=
%subplot(1,3,2);
SAP = [S.ShiftAddPeak];
Freq = [S.Freq_ipsi];
contourf(PH, Freq, fliplr(circshift(SAP ,Nh).'));
if makeYlabel, ylabel('Frequency (Hz)', 'fontsize',14); end
xlabel('IPD (cycles)', 'fontsize',14);
caxis(max(SAP(:))+[2*(Wmin-Wmax) 0]);
xlim([-0.5 0.5])
set(gca, 'xtick',-0.5:0.25:0.5, 'fontsize', 13, 'ytick', 0:250:1e4);

function local_powerfun(EarStr, pr);
PL_I = JLpowerlaw(pr.Uidx_I);
PL_C = JLpowerlaw(pr.Uidx_C);
PL_B = JLpowerlaw(pr.Uidx_B);
figure;
set(gcf,'units', 'normalized', 'position', [0.0906 0.613 0.313 0.244]+0.05*[randn(1,2) 0 0]);
DCfield = 'MeanDriv';
%DCfield = 'MeanSpont1';
plot(PL_I.Vinst-PL_I.(DCfield), PL_I.cond_spike_prob, 'b');
xplot(PL_C.Vinst-PL_C.(DCfield), PL_C.cond_spike_prob, 'r');
xplot(PL_B.Vinst-PL_B.(DCfield), PL_B.cond_spike_prob, 'k');
title(JLtitle(pr.Uidx_B), 'interpreter', 'none');
%xplot(AD.Vinst-CS.MeanDriv, AD.cond_spike_prob, 'b-<', 'markersize',4);

function [DIPratio, ha] = local_binvardip(Uidx, PanelChar, Xpos, PhaseShift);
PhaseShift= arginDefaults('PhaseShift',[0 0]);
T = JLanova2(Uidx);
S = JLdatastruct(Uidx);
CS = JLcycleStats(Uidx);
Phase = linspace(0,1-1/T.NsamPerCycle,T.NsamPerCycle);
PercBelow = 100*sum(T.var2d(:)<CS.VarSpont1)/numel(T.var2d)
SamShift = round(T.NsamPerCycle*PhaseShift([2 1]));
T.mean2d = circshift(T.mean2d, SamShift);
T.var2d = circshift(T.var2d, SamShift);
% =mean=
ha(1) = axes('position', [Xpos  0.65    1.1*[0.15    0.22]], 'fontsize', 12);
CC = contourf(Phase, Phase, T.mean2d);
% =var=
ha(2) = axes('position', [Xpos  0.25    1.1*[0.15    0.22]], 'fontsize', 12);
CC = contourf(Phase, Phase, T.var2d);
hold on
contour(Phase, Phase, T.var2d, CS.VarSpont1, 'color', 'w', 'linewidth', 3)
%contour(Phase, Phase, T.var2d, CS.VarSpont1, 'color', 'w', 'linewidth', 1)
%title('Across-beat variance', 'fontangle', 'italic');
local_mark(ha(1), PanelChar(1), -0.2, 0.3, 'backgroundcolor', 'none');
local_mark(ha(2), PanelChar(2), -0.2, 0.3, 'backgroundcolor', 'none');
DIPratio = min(T.var2d(:))/CS.VarSpont1;

function MM = local_featherplot(dt, Mn, St, CLR, lineCLR);
if numel(dt)==1, dt = [dt 0]; end;
lineCLR = arginDefaults('lineCLR', CLR);
Tm = Xaxis(Mn, dt(1))+dt(2);
Pstd = [Mn(:)+St(:) ; flipud(Mn(:)-St(:))];
Ptime = [Tm(:) ; flipud(Tm(:))];
gCLR = 0.85*(CLR+1)/max(CLR+1);
patch(Ptime, Pstd, gCLR, 'FaceAlpha', 0.5, 'edgecolor', gCLR);
xdplot(dt, Mn, 'linewidth', 2, 'color', lineCLR);
[MM(1), MM(2)] = minmax(Pstd);



