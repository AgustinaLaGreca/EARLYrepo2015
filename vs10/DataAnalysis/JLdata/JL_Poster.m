function [Yi, Yc, Yb] = JL_Poster(ifig, doPrint, Yi, Yc, Yb);
% JL_Poster(ifig, doPrint);
%
%     % =======panel 3: Periodicity analysis; MONO======
%     case 3.1, % raw trace of monaural tone
%     case 3.2, % zoom of raw trace of monaural tone
%     case 3.3, % cycle-averaged input waveform
%     case 3.4 % cycle-av : multiple peaks IPSI
%     case 3.5 % cycle-av : multiple peaks CONTRA
%     case 3.6 % cycle-av : multiple peaks CONTRA
%     case 3.7 % cycle-av : multiple peaks CONTRA
%     case 3.8 % multipeak: spike pred ~ 3.6
%     case 4.1 % binaural interaction: mean & var of I,C,B
%     case 4.2 % binaural interaction: I+C versus B
%     case 4.3 % elongated beat plot, upper half
%     case 4.4 % elongated beat plot, lower half
%     case 4.5 % I+C vs B
%     case 4.6 % bin spike pred
%     case 4.7 % bin spike pred MLTp

if nargin<2, doPrint=0; end

global Jb214_6  Jb214_4

if isempty(Jb214_6),
    Jb214_6 = local_read('RG10214', 6);
end
if isempty(Jb214_4),
    Jb214_4 = local_read('RG10214', 4);
end

StimMarker = struct('color', [1 0.6 0.6], 'linewidth', 4);
GrayText = struct('units', 'normal', 'color', 0.7*[1 1 1]);
APline = struct('color', [0 0.7 0.2], 'linewidth', 3);
APmarker = struct('marker', '*', 'color', [0 0.7 0.2], 'linewidth', 2, 'markersize', 10);
ThickGray = struct('color', 0.6*[1 1 1], 'linewidth', 2);
StdAx = @(h, varargin)set(h, 'fontsize', 12, 'Layer', 'top', varargin{:});
PU = CollectInStruct(APmarker, APline, GrayText, Jb214_6, StdAx, StimMarker, ThickGray);
switch ifig,
    % =======panel 3: Periodicity analysis; MONO======
    case 2.1  % spont trace
        D = readTKABF(Jb214_6.I50(7)); 
        set(gcf,'units', 'normalized', 'position', [-0.00234 0.585 1 0.323]);
        Toffset = 8800; % ms plot ofset
        rec = D.AD(1,1).samples;
        TT = timeaxis(rec, D.dt_ms, -Toffset);
        plot(TT, rec, 'k', 'linewidth', 1);
        xplot(TT, rec+pmask(betwixt(TT,[17.78 18.35])), APline);
        xplot(17.91, 6, APmarker);
        xlim([16 49.8]);
        %xlim([8977.8    8992.2]-Toffset);
        ylim([2.5 6.5]);
        %fenceplot(Jb214_6.I50(7).SPTraw-Toffset, ylim, 'r');
        text(0.75, 0.9, 'RG10214-6 / 6-13 (7)', GrayText, 'units', 'normaliz');
        StdAx(gca);
        xlabel('Time (ms)', 'fontsize', 13);
        ylabel('Potential (mV)', 'fontsize', 13);
        TracePlotInterface(gcf);
    case 3.1, % raw trace of monaural tone
        set(gcf,'units', 'normalized', 'position', [-0.00234 0.585 1 0.323]);
        Jb = Jb214_6;
        D = readTKABF(Jb.I50(4)); 
        rec = D.AD(1,1).samples;
        TT = timeaxis(rec,D.dt_ms);
        dplot(D.dt_ms, rec, 'k');
        xplot([500 6500], 2.3*[1 1], StimMarker); % stimulus
        xlim([0 8700])
        title('');
        text(0.75, 0.9, 'RG10214-6 / 300 Hz / 50 dB ipsi', GrayText);
        StdAx(gca);
        xlabel('Time (ms)', 'fontsize', 13);
        ylabel('Potential (mV)', 'fontsize', 13);
    case 3.2, % zoom of raw trace of monaural tone
        Jb = Jb214_6;
        set(gcf,'units', 'normalized', 'position', [-0.00234 0.585 1 0.323]);
        Jb = Jb214_6;
        D = readTKABF(Jb.I50(4)); 
        rec = D.AD(1,1).samples;
        TT = timeaxis(rec,D.dt_ms);
        dplot(D.dt_ms, rec, 'k', 'linewidth', 2);
        % highlight spike
        D = readTKABF(Jb.I50(4)); 
        rec = D.AD(1,1).samples;
        TT = timeaxis(rec,D.dt_ms);
        xlim([511.2000  530.8830])
        xplot(TT,rec+pmask(betwixt(TT,521.15,522.3)), APline);
        xplot(521.2902, 6.2729, '*', APmarker)
        xplot(TT, 2.3+0.7*cos(2*pi*TT*0.3), StimMarker); % stimulus
        ylim([1.2 6.8]);
        title('');
        text(0.75, 0.9, 'RG10214-6 / 300 Hz / 50 dB ipsi', GrayText);
        StdAx(gca);
        xlabel('Time (ms)', 'fontsize', 13);
        ylabel('Potential (mV)', 'fontsize', 13);
    case 3.3, % cycle-averaged input waveform
        Y = JLbeatVar(Jb214_6.I50(4)); % 400 Hz
        local_cyc_av_mon(Y, 1, PU);
    case 3.4 % cycle-av : multiple peaks IPSI
        Y = JLbeatVar(Jb214_6.I80(2)); % 200 Hz
        local_cyc_av_mon(Y, 1, PU);
    case 3.5 % cycle-av : multiple peaks CONTRA
        Y = JLbeatVar(Jb214_6.C80(2)); % 200 Hz
        local_cyc_av_mon(Y, 0, PU);
    case 3.6 % cycle-av : multiple peaks CONTRA
        Y = JLbeatVar(Jb214_4.I60(3)); % 300 Hz
        local_cyc_av_mon(Y, 1, PU, 0.45); % 0.4 ms circ shift
    case 3.7 % cycle-av : multiple peaks CONTRA
        Y = JLbeatVar(Jb214_4.C60(3)); % 300 Hz
        local_cyc_av_mon(Y, 0, PU, 0.45);
    case 3.8 % multipeak: spike pred ~ 3.6
        Y = JLbeatVar(Jb214_4.I60(3)); 
        P=JLpredictSpikes(Y, nan);
        Nbin = 25; Tshift = 0.45; % ms shift 
        BinCenter = ((1:Nbin)-0.5)/Nbin;
        PHipsi = P.SPTipsi/P.IpsiCycle;
        PHipsi = rem(PHipsi+Tshift/P.IpsiCycle,1);
        set(figure,'units', 'normalized', 'position', [0.0563 0.232 0.379 0.584]);
        hist(PHipsi, BinCenter);
        PredRate = P.Nspike*P.Ipsi_rate/P.dt_ipsi/P.IpsiCycle;
        PredRate = circshift(PredRate, round(Tshift/P.dt_ipsi));
        hist(PHipsi, BinCenter);
        local_colorpatch(gca, [0 0.1 0.8]);
        xdplot(P.dt_ipsi/P.IpsiCycle, PredRate, 'color', 0.6*[0 1 0], 'linewidth', 3);
        StdAx(gca);
        xlabel('Phase (cycle)', 'fontsize', 13);
        ylabel('Spike count', 'fontsize', 13);
        text(0.43, 0.9, 'RG10214-4 / 300 Hz / 60 dB ipsi', GrayText);
    case 4.1 % binaural interaction: mean & var of I,C,B
        local_binint(Jb214_6.I50(4), Jb214_6.C50(4), Jb214_6.B50(4), PU);
    case 4.2 % binaural interaction: I+C versus B
        [Yi, Yc, Yb] = local_binint(Jb214_6.I50(4), Jb214_6.C50(4), Jb214_6.B50(4), 'noplot');
        fh2 = figure;
        set(fh2,'units', 'normalized', 'position', [0.532 0.166 0.448 0.529]);
        plot(Yi.ipsi_MeanWave, Yb.ipsi_MeanWave, 'b', 'linewidth', 2);
        xplot(Yc.contra_MeanWave, Yb.contra_MeanWave, 'color', [0.8 0 0], 'linewidth', 2);
        StdAx(gca);
        xlim([3 5]); 
        ylim([2.5 4.5]);
        set(gca, 'xtick', -10:0.5:10)
        set(gca, 'ytick', -10:0.5:10)
        axis square; 
        grid on
        xlabel('Monaural response (mV)', 'fontsize', 13);
        ylabel('Contribution to binaural response (mV)', 'fontsize', 13);
        text(4.35, 3.85, 'ipsi', 'fontweight', 'bold', 'fontsize', 12, 'color', 'b');
        text(3.7, 3.7, 'contra', 'fontweight', 'bold', 'fontsize', 12, 'color', [0.8 0 0]);
    case 4.3 % elongated beat plot, upper half
        if nargin<3,
            [Yi, Yc, Yb] = local_binint(Jb214_6.I50(4), Jb214_6.C50(4), Jb214_6.B50(4), 'noplot');
        end
        [TTb, Ri, Rc, Rb] = local_interp(Yi,Yc,Yb);
        set(gcf,'units', 'normalized', 'position', [0 0.799 1 0.0725])
        %dsize(TTi, TTc, TTb, Ri, Rc, Rb)
        plot(TTb, Ri, 'b'); 
        xplot(TTb, Rc, 'r'); 
        Yoffset = 1.15;
        xplot(TTb, Rb+Yoffset, 'k', 'linewidth', 1.5);
        xplot(TTb, Ri+Rc+Yoffset, 'm');
        set(gca,'xtick', [], 'ytick', []);
        ylim([-0.5 2.3]);
        xlim([0 125])
    case 4.4 % elongated beat plot, lower half
        if nargin<3,
            [Yi, Yc, Yb] = local_binint(Jb214_6.I50(4), Jb214_6.C50(4), Jb214_6.B50(4), 'noplot');
        end
        [TTb, Ri, Rc, Rb] = local_interp(Yi,Yc,Yb);
        set(gcf,'units', 'normalized', 'position', [0 0.799 1 0.0725])
        %dsize(TTi, TTc, TTb, Ri, Rc, Rb)
        plot(TTb, Ri, 'b'); 
        xplot(TTb, Rc, 'r'); 
        Yoffset = 1.15;
        xplot(TTb, Rb+Yoffset, 'k', 'linewidth', 1.5);
        xplot(TTb, Ri+Rc+Yoffset, 'm');
        set(gca,'xtick', [], 'ytick', []);
        ylim([-0.5 2.3]);
        xlim([125 250])
    case 4.5 % I+C vs B
        if nargin<3,
            [Yi, Yc, Yb] = local_binint(Jb214_6.I50(4), Jb214_6.C50(4), Jb214_6.B50(4), 'noplot');
        end
        [TTb, Ri, Rc, Rb] = local_interp(Yi,Yc,Yb);
        set(gcf,'units', 'normalized', 'position', [0.28 0.226 0.502 0.645]);
        plot(Ri+Rc, Rb, 'k', 'linewidth', 1);
        StdAx(gca);
        xlabel('Prediction I+C (mV)', 'fontsize', 13);
        ylabel('Measured binaural response (mV)', 'fontsize', 13);
        xlim([-0.8 1.1]); ylim(xlim);
        set(gca,'xtick', -1:0.5:2,'ytick', -1:0.5:2);
        axis square;
        grid on;
    case 4.6 % bin spike pred
        if nargin<3,
            [Yi, Yc, Yb] = local_binint(Jb214_6.I50(4), Jb214_6.C50(4), Jb214_6.B50(4), 'noplot');
        end
        Tshift = 0;
        P=JLpredictSpikes(Yb, nan); 
        Nbin = 25; Tshift = 0.45; % ms shift 
        BinCenter = ((1:Nbin)-0.5)/Nbin;
        PHipsi = P.SPTipsi/P.IpsiCycle;
        PHipsi = rem(PHipsi+Tshift/P.IpsiCycle,1);
        PHcontra = P.SPTcontra/P.ContraCycle;
        PHcontra = rem(PHcontra+Tshift/P.ContraCycle,1);
        PHbin = P.SPTbin/P.BeatCycle;
        PHbin = rem(PHbin+(Tshift+125)/P.BeatCycle,1);
        set(gcf,'units', 'normalized', 'position', [0.0203 0.584 0.932 0.26]);
        subplot(1,3,1);
        PredRate = P.Nspike*P.Ipsi_rate/P.dt_ipsi/P.IpsiCycle;
        PredRate = 0.7*circshift(PredRate, round(Tshift/P.dt_ipsi));
        hist(PHipsi, BinCenter);
        local_colorpatch(gca, [0 0.1 0.8]);
        xdplot(P.dt_ipsi/P.IpsiCycle, PredRate, 'color', 0.6*[0 1 0], 'linewidth', 3);
        StdAx(gca); title('IPSI 400 Hz');
        xlabel('Phase (cycle)', 'fontsize', 13);
        ylabel('Spike count', 'fontsize', 13);
        subplot(1,3,3);
        PredRate = P.Nspike*P.Contra_rate/P.dt_contra/P.ContraCycle;
        PredRate = 1.25*circshift(PredRate, round(Tshift/P.dt_contra));
        hist(PHcontra, BinCenter);
        local_colorpatch(gca, [0.9 0.05 0]);
        xdplot(P.dt_contra/P.ContraCycle, PredRate, 'color', 0.6*[0 1 0], 'linewidth', 3);
        StdAx(gca); title('CONTRA 404 Hz');
        xlabel('Phase (cycle)', 'fontsize', 13);
        subplot(1,3,2); %===================================
        PredRate = P.Nspike*P.APrate/P.dt_bin/P.BeatCycle;
        PredRate = 150*circshift(PredRate, round((Tshift+125)/P.dt_bin));
        PredRate = smoothen(PredRate, 15, P.dt_bin);
        PredRate = PredRate - min(PredRate);
        hist(PHbin, BinCenter);
        xdplot(P.dt_bin/P.BeatCycle, PredRate, 'color', 0.6*[0 1 0], 'linewidth', 3);
        local_colorpatch(gca, [0.8 0.0 0.9]);
        StdAx(gca); title('BIN BEAT 4 Hz');
        xlabel('Phase (cycle)', 'fontsize', 13);
        %ylabel('Spike count', 'fontsize', 13);
        ylim([0 35]);
        text(0.1, 0.9, 'RG10214-6 / 400 Hz / 50 dB ipsi', GrayText);
    case 4.7 % bin spike pred MLTp
        Yb = JLbeatVar(Jb214_6.B80(2),0); 
        Tshift = 0;
        P=JLpredictSpikes(Yb, nan); 
        Nbin = 25; Tshift = 0.45; % ms shift 
        BinCenter = ((1:Nbin)-0.5)/Nbin;
        PHipsi = P.SPTipsi/P.IpsiCycle;
        PHipsi = rem(PHipsi+Tshift/P.IpsiCycle,1);
        PHcontra = P.SPTcontra/P.ContraCycle;
        PHcontra = rem(PHcontra+Tshift/P.ContraCycle,1);
        PHbin = P.SPTbin/P.BeatCycle;
        PHbin = rem(PHbin+(Tshift+125)/P.BeatCycle,1);
        set(gcf,'units', 'normalized', 'position', [0.0203 0.584 0.932 0.26]);
        subplot(1,3,1);
        PredRate = P.Nspike*P.Ipsi_rate/P.dt_ipsi/P.IpsiCycle;
        PredRate = 1.3*circshift(PredRate, round(Tshift/P.dt_ipsi));
        hist(PHipsi, BinCenter);
        local_colorpatch(gca, [0 0.1 0.8]);
        xdplot(P.dt_ipsi/P.IpsiCycle, PredRate, 'color', 0.6*[0 1 0], 'linewidth', 3);
        StdAx(gca); title('IPSI 200 Hz');
        xlabel('Phase (cycle)', 'fontsize', 13);
        ylabel('Spike count', 'fontsize', 13);
        subplot(1,3,3);
        PredRate = P.Nspike*P.Contra_rate/P.dt_contra/P.ContraCycle;
        PredRate = 1.9*circshift(PredRate, round(Tshift/P.dt_contra));
        hist(PHcontra, BinCenter);
        local_colorpatch(gca, [0.9 0.05 0]);
        xdplot(P.dt_contra/P.ContraCycle, PredRate, 'color', 0.6*[0 1 0], 'linewidth', 3);
        StdAx(gca); title('CONTRA 204 Hz');
        xlabel('Phase (cycle)', 'fontsize', 13);
        subplot(1,3,2); %===================================
        PredRate = P.Nspike*P.APrate/P.dt_bin/P.BeatCycle;
        PredRate = 160*circshift(PredRate, round((Tshift+125)/P.dt_bin));
        PredRate = smoothen(PredRate, 10, P.dt_bin);
        PredRate = PredRate - min(PredRate);
        hist(PHbin, BinCenter);
        xdplot(P.dt_bin/P.BeatCycle, PredRate, 'color', 0.6*[0 1 0], 'linewidth', 3);
        local_colorpatch(gca, [0.8 0.0 0.9]);
        StdAx(gca); title('BIN BEAT 4 Hz');
        xlabel('Phase (cycle)', 'fontsize', 13);
        %ylabel('Spike count', 'fontsize', 13);
        ylim([0 35]);
        text(0.1, 0.9, 'RG10214-6 / 200 Hz / 80 dB ipsi', GrayText);
    otherwise,
        error('unkown panel.fig');
end

if doPrint,
    switch upper(CompuName),
        case 'CEL',
            PosterDir = 'C:\D_Drive\doc\Meetings and Talks\ARO\ARO-2011\PosterFigs';
        case 'SIUT',
            PosterDir = '\\Ee1285a\ARO-2011$\PosterFigs';
    end
    set(gcf, 'paperpositionmode', 'auto');
    PanelStr = num2str(floor(ifig));
    FigStr = num2str(round(10*(rem(ifig,1))));
    FN = fullfile(PosterDir, ['panel_' PanelStr '_fig_' FigStr ])
    print(gcf, '-dpng', '-r600', FN);
end

%===========================================================

function Jb = local_read(ExpID, icell);
if isequal({'RG10214', 6}, {ExpID, icell}),
    JLreadBeats(ExpID, icell);   % 200 Hz & up, 100-Hz steps
    Jb = struct('B80', JbB_80dB, 'I80', JbI_80dB, 'C80', JbC_80dB, ...
        'B50', JbB_50dB, 'I50', JbI_50dB, 'C50', JbC_50dB);
elseif isequal({'RG10214', 4}, {ExpID, icell}),
    JLreadBeats(ExpID, icell);
    Jb = struct('B60', JbB_60dB, 'I60', JbI_60dB, 'C60', JbC_60dB);
else,
    error('don''t know what to read')
end

function local_featherplot(dt, Mn, Vr, CLR, lineCLR);
lineCLR = arginDefaults('lineCLR', CLR);
Tm = Xaxis(Mn, dt);
St = sqrt(Vr);
Pstd = [Mn(:)+St(:) ; flipud(Mn(:)-St(:))];
Ptime = [Tm(:) ; flipud(Tm(:))];
gCLR = 0.85*(CLR+1)/max(CLR+1);
set(gcf,'renderer','zbuffer');
patch(Ptime, Pstd, gCLR, 'FaceAlpha', 0.5, 'edgecolor', gCLR);
xdplot(dt, Mn, 'linewidth', 2, 'color', lineCLR);

function local_cyc_av_mon(Y, Ipsi, PU, tshift);
if nargin<4, tshift=0; end
Nshift = round(tshift./[Y.ipsi_dt Y.contra_dt]);
[Y.ipsi_MeanWave, Y.ipsi_VarWave] = ...
    deal(circshift(Y.ipsi_MeanWave, Nshift(1)), circshift(Y.ipsi_VarWave, Nshift(1)));
[Y.contra_MeanWave, Y.contra_VarWave] ...
    = deal(circshift(Y.contra_MeanWave, Nshift(2)), circshift(Y.contra_VarWave, Nshift(2)));
set(figure,'units', 'normalized', 'position', [0.163 0.0887 0.401 0.801]);
subplot(2,1,2); % mean
set(gca,'layer', 'top')
if Ipsi,
     local_featherplot(Y.ipsi_dt, Y.ipsi_MeanWave, Y.ipsi_VarWave, [0 0 1]);
else,
    local_featherplot(Y.contra_dt, Y.contra_MeanWave, Y.contra_VarWave, [1 0 0]);
end
%local_featherplot(Y.contra_dt, Y.contra_MeanWave, Y.contra_VarWave, [1 0 0]);
xplot(xlim, [1 1]*Y.SpontMean, PU.ThickGray);
ylabel('Cycle-avg. potential (mV)', 'fontsize', 13);
xlabel('Time (ms)', 'fontsize', 13);
PU.StdAx(gca, 'box', 'off');
set(gca,'layer', 'top');
xplot(xlim, min(ylim)*[1 1], 'k', 'linewidth',1);

subplot(2,1,1); % var

if Ipsi,
    dplot(Y.ipsi_dt, Y.ipsi_VarWave, 'b', 'linewidth', 2);
    xplot(xlim, [1 1]*Y.SpontVar, PU.ThickGray)
    %xdplot(Y.contra_dt, Y.contra_VarWave, 'r');
    xplot(xlim, mean(Y.ipsi_VarWave)*[1 1], 'b--', 'linewidth', 2)
    xplot(xlim, mean(Y.contra_VarWave)*[1 1], 'r--', 'linewidth', 2)
    title('IPSI', 'fontsize', 14);
else,
    dplot(Y.contra_dt, Y.contra_VarWave, 'r', 'linewidth', 2);
    xplot(xlim, [1 1]*Y.SpontVar, PU.ThickGray)
    xdplot(Y.ipsi_dt, Y.ipsi_VarWave, 'b');
    xplot(xlim, mean(Y.contra_VarWave)*[1 1], 'r--', 'linewidth', 2)
    xplot(xlim, mean(Y.ipsi_VarWave)*[1 1], 'b:')
    title('CONTRA', 'fontsize', 14);
end
PU.StdAx(gca, 'box', 'off');
ylim([0 max(ylim)]);
ylabel('Across-cycle variance (mV^2)', 'fontsize', 13);
set(gca,'layer', 'top');
xplot(xlim, min(ylim)*[1 1], 'k', 'linewidth',1);


function local_colorpatch(h, CLR);
hp = findobj(h, 'type', 'patch');
set(hp, 'facecolor', CLR);

% fh1 = figure;
% set(fh1,'units', 'normalized', 'position', [0.0391 0.471 0.658 0.456])
% figname(fh1, [mfilename ' mean&var vs time'])

function [Yi, Yc, Yb] = local_binint(Yi, Yc, Yb, PU);
if ~isfield(Yi, 'bin_MeanWave'), % beatplot format; convert to beatVar format
    Yi = JLbeatVar(Yi);
    Yc = JLbeatVar(Yc);
    Yb = JLbeatVar(Yb);
end
if isequal('noplot', PU), return; end
set(gcf,'units', 'normalized', 'position', [0.28 0.185 0.652 0.686])
[Rmax_ipsi, tau_ipsi]    = local_corr(Yb.ipsi_dt,   Yb.ipsi_MeanWave, Yi.ipsi_MeanWave);
[Rmax_contra, tau_contra]= local_corr(Yb.contra_dt, Yb.contra_MeanWave, Yc.contra_MeanWave);
Istr = sprintf('\\rho=%1.4f; delay %i \\mus', Rmax_ipsi, tau_ipsi);
Cstr = sprintf('\\rho=%1.4f, delay %i \\mus', Rmax_contra, tau_contra);
% var
h1 = subplot(2,3,1); PU.StdAx(h1); 
dplot(Yi.ipsi_dt, Yi.ipsi_VarWave, 'b', 'linewidth',2);
xlim([0 Yb.Period1]);
xplot(xlim, [1 1]*Yi.SpontVar, PU.ThickGray);
ylabel('Across-cycle variance (mV^2)', 'fontsize', 13);
title(sprintf('IPSI %i Hz', round(Yi.Freq1)));
h2 = subplot(2,3,2); PU.StdAx(h2);
dplot(Yb.ipsi_dt, Yb.ipsi_VarWave, 'b', 'linewidth' ,2);
xdplot(Yb.contra_dt, Yb.contra_VarWave, 'r', 'linewidth',2);
xplot(xlim, [1 1]*Yb.SpontVar, PU.ThickGray);
title(sprintf('BIN. BEAT %i Hz', round(abs(Yb.Freq2-Yb.Freq1))));
h3 = subplot(2,3,3); PU.StdAx(h3);
dplot(Yc.ipsi_dt, Yc.contra_VarWave, 'r', 'linewidth' ,2);
xplot(xlim, [1 1]*Yc.SpontVar, PU.ThickGray);
ylim(h1, [0 max([ylim(h1) ylim(h2) ylim(h3)])]);
linkaxes([h1 h2 h3], 'xy')
set([h2 h3], 'yticklabel',{});
title(sprintf('CONTRA %i Hz', round(Yc.Freq2)));
% mean
h4 = subplot(2,3,4); PU.StdAx(h4);
dplot(Yi.ipsi_dt, Yi.ipsi_MeanWave, 'b', 'linewidth' ,2);
xlim([0 Yb.Period1]);
%xplot(xlim, [1 1]*Yi.SpontMean, ThickGray);
xlabel('Time (ms)', 'fontsize', 13);
ylabel('Cycle-averaged potential (mV)', 'fontsize', 13);
h5 = subplot(2,3,5); PU.StdAx(h5);
dplot(Yb.ipsi_dt, Yb.ipsi_MeanWave, 'b', 'linewidth' ,2);
xdplot(Yb.contra_dt, Yb.contra_MeanWave, 'r', 'linewidth' ,2);
text(0.16, 0.9, Istr, 'units', 'norm', 'color', 'b', 'fontweight', 'bold');
text(0.16, 0.8, Cstr, 'units', 'norm', 'color', 'r', 'fontweight', 'bold');
%xplot(xlim, [1 1]*Yb.SpontMean, ThickGray);
xlabel('Time (ms)', 'fontsize', 13);
h6 = subplot(2,3,6); PU.StdAx(h6);
dplot(Yc.ipsi_dt, Yc.contra_MeanWave, 'r', 'linewidth' ,2);
%xplot(xlim, [1 1]*Yc.SpontMean, ThickGray);
YL = [min([ylim(h4) ylim(h5) ylim(h6)]) max([ylim(h4) ylim(h5) ylim(h6)])];
ylim(h4, [YL(1) mean(YL)+0.8*diff(YL)]);
linkaxes([h4 h5 h6], 'xy')
% legend(h5,{sprintf('R=%1.3f; delay %i \\mus', Rmax_ipsi, tau_ipsi), ...
%     sprintf('R=%1.3f, delay %i \\mus', Rmax_contra, tau_contra)}, ...
%     'location', 'North', 'fontsize', 9);
xlabel('Time (ms)', 'fontsize', 13);


function [rhomax, tau]=local_corr(dt, BinW, MonW);
Nmax = numel(BinW);
[rhomax, itau] = maxcorr(repmat(BinW, 10,1), repmat(MonW, 10,1), Nmax);
tau = round(1e3*dt*itau); % in mus

function [TTb, Ri, Rc, Rb] = local_interp(Yi,Yc,Yb);
Ri = repmat(Yi.ipsi_MeanWave-mean(Yi.ipsi_MeanWave), Yi.Freq1/4,1);
Rc = repmat(Yc.contra_MeanWave-mean(Yc.contra_MeanWave), Yc.Freq2/4,1);
Rb = Yb.bin_MeanWave-mean(Yb.bin_MeanWave);
TTi = Xaxis(Ri, Yi.ipsi_dt);
TTc = Xaxis(Rc, Yc.contra_dt);
TTb = Xaxis(Rb, Yb.bin_dt);
Ri = interp1(TTi, Ri, TTb);
Rc = interp1(TTc, Rc, TTb);











