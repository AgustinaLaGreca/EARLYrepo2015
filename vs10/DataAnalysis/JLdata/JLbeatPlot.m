function H = JLbeatPlot(S, iselect, flag);
% JLbeatPlot - plot binaural beat analysis analysis from JLbeat
%   JLbeatPlot(S) plots the info in S, the return variable of JLbeat, and
%   returns the plot handles of the axes.
%
%   JLbeatPlot(S, K, Flag) only shows the type K plot. 
%     1: recording + stimulus waveforms
%     2: spectrum
%     3: periodic analysis of recording 
%     4: I/O curve between linear model & recording 
%     5: spike-triggered average of truncated waveforms
%     6: spike histograms
%     7: AP snippets
%   The Flag string realizes additional plot settings:
%     'x3': xlimit([0 30]) in plot 3.
%
%
%   JLbeatPlot({S1 S2 ..}, K) plot S1, S2, ... in different figures.
%   If there are three Ss, figs are aligned on screen for a  
%   Bin/Ipsi/Contra comparison.
%
%   JLbeatPlot([S1 S2 ..], K) pools S1, S2, .. and plots the pooled data.
%
%   See also JLbeat.


if nargin<2, iselect=0; end % default: all plot types
if nargin<3, flag=''; end
Xtra = localXtra(flag);

if isequal(0,iselect),
    iselect = 1:7;
end

if isstruct(S) && numel(S)>1, % pool the data
    S = JLbeatPool(S);
end
if iscell(S) && numel(S)>1, % recursive
    % first visit all elements of S to see if pooling is needed
    for ii=1:numel(S),
        if numel(S{ii})>1,
            S{ii} = JLbeatPool(S{ii});
        end
    end
    % convert cell array of struct -> struct array
    for ii=2:numel(S),
        S{1} = [S{1} S{ii}];
    end
    S = S{1};
    H = [];
    [S.Yoffset] = deal(max([S.Yoffset]));
    for iplot=1:numel(iselect),
        for idata=1:numel(S),
            h = JLbeatPlot(S(idata), iselect(iplot), flag);
            H = [H h];
        end
        linkaxes(H,'xy');
    end
    if isequal(3, numel(S)) && numel(iselect)==1,
        % probably a bin/ipsi/contra comparison. Arrange figures.
        fh = []; for hh=H(:).', fh = [fh get(hh,'parent')]; end; fh = unique(fh);
        switch iselect,
            case {1 3}, 
                BasePos = [0.0169 0.0375 0.943 0.335];
                PosShift = [0 0.3 0 0];
            case 2, 
                BasePos = [0.261 0.07 0.683 0.267];
                PosShift = [0 0.3 0 0];
            case {4 5 7}, 
                BasePos = [0.559 0.07 0.315 0.266];
                PosShift = [0 0.3 0 0];
            case {6}, 
                BasePos = [0.38 0.03 0.149 0.308];
                PosShift = [0 0.3 0 0];
        end
        set(fh(3),'units', 'normalized', 'position', BasePos + 2*PosShift); 
        set(fh(2),'units', 'normalized', 'position', BasePos + PosShift);
        set(fh(1),'units', 'normalized', 'position', BasePos);
    end
    if nargout<1, clear('H'); end % only return handles if explicitly requested
    return;
end
%====single S from here=========
S = JLgetBeatStruct(S);
H = [];

if ismember(1, iselect),
    % plot recording together with waveforms
    fh1 = figure;
    set(fh1,'units', 'normalized', 'position', [0.0188 0.462 0.98 0.438]);
    dplot(S.dt, S.R0, 'k', 'linewidth', 2);
    xdplot(S.dt, -S.Yoffset + S.Stim1, 'b');
    xdplot(S.dt, -S.Yoffset + S.Stim2, 'r');
    ylim(ylim);
    xlabel('Time (ms)');
    title(S.TTT);
    legend({'Recording' 'Left-ear stimulus' 'Left-ear stimulus'});
    TracePlotInterface(fh1);
    ax1 = get(fh1, 'CurrentAxes');
    H = [H ax1];
    JLfigmenu(fh1,S);
end

if ismember(2, iselect),
    % spectral plots
    fh2 = figure;
    set(fh2,'units', 'normalized', 'position', [0.324 0.564 0.491 0.35]);
    FreqAx = timeaxis(S.R0,S.df);
    Spec = A2dB(abs(fft(S.R0)));
    Spec(1) = nan; % ignore DC
    plot(FreqAx, Spec, 'k');
    xlim([0 3000]);
    xlabel('Frequency (Hz)');
    for icomp=numel(S.Sb):-1:1, % in order of descending order ;)
        sb = S.Sb(icomp);
        iord = sb.order;
        isHarm = isequal(sb.order, max(sb.weight)); % harmonic 
        idx = 1 + round([sb.freq]/S.df);
        hdp(iord)=xplot(FreqAx(idx), Spec(idx), '*', lico(iord));
        if isHarm, 
            xplot(FreqAx(idx), Spec(idx), 'o', lico(iord))
        end
        LegStr{iord} = ['order = ' num2str(iord)];
    end
    legend(hdp, LegStr);
    ylim(max(ylim)+[-80 0]);
    title(S.TTT);
    H = [H get(fh2, 'CurrentAxes')];
    JLfigmenu(fh2,S);
end


if ismember(3, iselect),
    % spectral plots
    fh3 = figure;
    % decomposition in ipsi & contra contributions
    ttime = (0:numel(S.Y1)-1).'*S.dt;
    % original recording and comb-filtered versions
    hR0 = dplot(S.dt, S.R0, 'k', 'linewidth', 2);
    hYc = xdplot(S.dt, S.Yc, 'g');
    hY3 = xdplot(S.dt, S.Y3, 'm');
    % decomposition into ipsi & contra
    h0 = xdplot(S.dt, 0*S.R0-S.Yoffset, 'k');
    xdplot(S.dt, S.Y1-S.Yoffset, 'b');
    xdplot(S.dt, S.Y2-S.Yoffset, 'r');
    pm1 = pmask(~rem(floor(ttime/S.Period1),2));
    pm2 = pmask(~rem(floor(ttime/S.Period2),2));
    hY1 = xdplot(S.dt, S.Y1+pm1-S.Yoffset,'b', 'linewidth',1.5);
    hY2 = xdplot(S.dt, S.Y2+pm2-S.Yoffset,'r', 'linewidth',1.5);
    % residues
    xdplot(S.dt, 0*S.R0+S.Yoffset, 'k');
    xdplot(S.dt, (S.R0-S.Yc)+S.Yoffset, 'g');
    xdplot(S.dt, (S.R0-S.Y3)+S.Yoffset, 'm');
    % spike times
    fenceplot(S.SPT, 0.2*S.Yoffset*[-1 1], 'color', [0 0.7 0], 'linewidth', 2);
    % various
    set(gcf,'units', 'normalized', 'position', [0.0188 0.462 0.98 0.438]);
    if Xtra.xlim3, xlim([0 30]); end
    ylim(ylim); % fix Y limits so they won't jump when X-scrolling
    xlabel('Time (ms)');
    title([S.TTT ' \rho=' num2str(S.AC)]);
    legend([hR0 hY1 hY2 hYc hY3], {'Beat average', 'Left average', 'Right average', 'Sum prediction', 'Interaction prediction'})
    TracePlotInterface(fh3);
    if exist('fh1', 'var'),
        linkaxes([gca ax1], 'x');
    end
    H = [H get(fh3, 'CurrentAxes')];
    JLfigmenu(fh3,S);
end

if ismember(4, iselect),
    % X-Y plot of linear model vs data
    fh4 = figure;
    set(fh4,'units', 'normalized', 'position', [0.758 0.18 0.228 0.246]);
    plot(S.Yc, S.R0, '.');
    axis equal
    axis square
    grid on
    xplot(xlim, xlim, 'r-');
    H = [H get(fh4, 'CurrentAxes')];
    JLfigmenu(fh4,S);
end


if ismember(5, iselect),
    % spike trig averages
    fh5 = figure; 
    set(gcf,'units', 'normalized', 'position', [0.451 0.133 0.315 0.293]);
    TimeWin =  [-2 1];
    SPTa0 = SpikeTrigAv(S.dt, S.R0, S.SPT, TimeWin);
    SPTa1 = SpikeTrigAv(S.dt, S.Y1, S.SPT, TimeWin);
    SPTa2 = SpikeTrigAv(S.dt, S.Y2, S.SPT, TimeWin);
    Tw = timeaxis(SPTa0, S.dt, TimeWin(1));
    plot(Tw, SPTa0, 'k');
    xplot(Tw, SPTa1, 'b');
    xplot(Tw, SPTa2, 'r');
    title([num2str(numel(S.SPT)) ' spikes'])
    grid on
    H = [H get(fh5, 'CurrentAxes')];
    JLfigmenu(fh5,S);
end

if ismember(6, iselect),
    % period histograms
    fh6 = figure; 
    set(fh6,'units', 'normalized', 'position', [0.313 0.109 0.149 0.307])
    SPT1 = mod(S.SPT, S.Period1)/S.Period1;
    SPT2 = mod(S.SPT, S.Period2)/S.Period2;
    Nspike = numel(S.SPT);
    Tbeat = 1e3/S.Fbeat; % beat period in ms
    SPTb = mod(S.SPT, Tbeat)/Tbeat;
    Nspike = numel(SPT1);
    Nhist = min(50,round(Nspike/7));
    Nhist = max(Nhist, 10);
    Ph = ((1:Nhist)-0.5)/Nhist;
    R1 = abs(vectorstrength(S.SPT,S.Freq1));
    R2 = abs(vectorstrength(S.SPT,S.Freq2));
    Rb = abs(vectorstrength(S.SPT,S.Fbeat));
    ah1=subplot(3,1,1);
    hist(SPT1, Ph);
    text(0.5, 0.8, ['ipsi: R=' num2str(R1,2)], 'units', 'normalized');
    ah2=subplot(3,1,2);
    hist(SPT2, Ph);
    text(0.5, 0.8, ['contra: R=' num2str(R2,2)], 'units', 'normalized');
    ah3=subplot(3,1,3);
    hist(SPTb, Ph);
    text(0.5, 0.8, ['bin: R=' num2str(Rb,2)], 'units', 'normalized');
    xlabel('Phase (cycle)');
    H = [H ah1 ah2 ah3];
    disp(['% ' num2str(S.Freq1) ' Hz  R = [' num2str(R1,2) ' ' num2str(R2,2) ' ' num2str(Rb,2) '] ' ...
        num2str(Nspike),  ' spikes; rho=' num2str(deciRound(S.AC,2),3), ';'])
end

if ismember(7, iselect),
    % spike snippets
    fh7 = figure;
    set(fh7,'units', 'normalized', 'position', [0.0219 0.188 0.295 0.252]);
    if ~isempty(S.Snip), 
        % Iap = categorize(S.APslope, 5);
        % catplot(S.Tsnip, S.Snip, Iap);
        plot(S.Tsnip, S.Snip);
    end;
    xlabel('Time (ms)');
    ylabel('V_{rec} (mV)');
    title([num2str(size(S.Snip,2)) ' Spikes']);
    grid on;
    H = [H get(fh7, 'CurrentAxes')];
    JLfigmenu(fh7,S);
end


if exist('fh3','var'), figure(fh3); end


%===========================
function Xtra = localXtra(flag);
Xtra = struct('xlim3', 0);
F = Words2cell(flag,'/');
for ii=1:numel(F),
    switch lower(F{ii}),
        case 'x3';
            Xtra.xlim3 = 1;
    end
end



