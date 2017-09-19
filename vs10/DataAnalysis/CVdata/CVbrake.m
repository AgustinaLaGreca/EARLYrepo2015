function CVbrake(ExpName, ha);
% CVbrake analyze sudden slowdown of traveling wave
%    CVbrake(ExpName); 
%    CVbrake(ExpName, hax); 
ha = arginDefaults('ha',[]);

switch ExpName
    case 'all',
        figure;
        set(gcf,'units', 'normalized', 'position', [0.216 0.365 0.429 0.54], ...
            'PaperOrientation', 'landscape', 'paperpositionmod', 'auto');
        for ii=1:8, ah(ii) = subplot(2,4,ii); end
        CVbrake('RG12401', ah([1 5])); 
        CVbrake('RG12422', ah([2 6])); 
        CVbrake('RG12423', ah([3 7]));
        CVbrake('RG12424', ah([4 8]));
        for ii=[2:4 6:8], ylabel(ah(ii),''); end
        for ii=1:4, xlabel(ah(ii),''); end
        return;
    case 'RG12401',
        iSPL = [2 2 0];
        [FN1 FN2] = deal('L6', 'L9');
        Fwin = [5.5 13]; % kHz view window 
        cdelay = 75; % us
        L1 = [5.9219    0.376];
        L2 = [10.4219    0.624];
        H1 = [9.7031    0.66598];
        H2 = [12.9219    0.2332];
        YL = [0.2 0.7];
        Fx = 10.05;
    case 'RG12422',
        iSPL = [2 2 0];
        [FN1 FN2] = deal('L10', 'L16');
        Fwin = [10 18]; % kHz view window 
        cdelay = 75; % us
        L1 = [10.1523    0.6174];
        L2 = [14.5352    0.8596];
        H1 = [14.0195    0.8414];
        H2 = [17.9336    0.4430];
        YL = [0.4  0.9];
        Fx = 13.7;
    case 'RG12423',
        iSPL = [2 0 0];
        [FN1 FN2] = deal('L10', 'L16');
        Fwin = [10 17]; % kHz view window 
        cdelay = 95; % us
        L1 = [10.0638    0.8760];
        L2 = [13.6354    1.0865];
        H1 = [12.8854    1.0906];
        H2 = [16.5898    0.8052];
        YL = [0.75  1.12];
        Fx = 12.8;
    case 'RG12424',
        iSPL = [2 2 0];
        [FN1 FN2] = deal('L10', 'L16');
        Fwin = [11.5 18]; % kHz view window 
        cdelay = 125; % us
        L1 = [11.8    1.3309];
        L2 = [15.2    1.7141];
        H1 = [14.0645    1.6734];
        H2 = [17.8340    1.24];
        YL = [1.2 1.8];
        Fx = 14.05;
end
qq = twobeadstuff('getapple', ExpName);
abw = qq.Lw.A2B(iSPL(1));
if iSPL(2)>0,
    ab1 = qq.(FN1).A2B(iSPL(2));
else,
    ab1 = [];
end
if iSPL(3)>0,
    ab2 = qq.(FN2).A2B(iSPL(3));
else,
    ab2 = [];
end

if isempty(ha), % new figure with 2 subplots
    set(gcf,'units', 'normalized', 'position', [0.193 0.374 0.258 0.525]);
    ha(1) = subplot(2,1,1);
    ha(2) = subplot(2,1,2);
end
PLA = @(m,c)struct('marker', m, 'color',c,'markerfacecolor',c)
%====ampl==
axes(ha(1));
local_Gainplot(ab1, PLA('^', [0 0.8 0]));
local_Gainplot(ab2, PLA('v', [0.8 0 0]));
local_Gainplot(abw, PLA('o', [0 0 0.8]));
xlim(Fwin); grid on
title(ExpName);
ylim(ylim);
%====phase==
axes(ha(2));
local_Phaseplot(ab1, cdelay, PLA('^', [0 0.8 0]));
local_Phaseplot(ab2, cdelay, PLA('v', [0.8 0 0]));
local_Phaseplot(abw, cdelay, PLA('o', [0 0 0.8]));
xlim(Fwin); grid on
ylim(YL);
PG = local_trianplot(L1,L2,H1,H2, cdelay);

axes(ha(1));
text(min(xlim)+0.02*diff(xlim), min(ylim)+0.05*diff(ylim), ...
    ['Gain pred ' num2str(PG) ' dB'], 'fontweight', 'bold', 'backgroundcolor', 'w')

%local_markGain(ha(1), abw, Fx, PG);
% ===========================
function local_Gainplot(ab, plotArgs);
if isempty(ab), return; end
xplot(ab.Fprim/1e3, ab.Gain+pmask(ab.Alpha<=0.001), plotArgs, 'markersize', 5, 'linestyle', '-');
ylabel('Gain (dB)')
set(gca,'xtick',0:2:30);

function local_Phaseplot(ab, cdelay, plotArgs);
if isempty(ab), return; end
freq = ab.Fprim/1e3; % kHz
PH = ab.Phase+pmask(ab.Alpha<=0.001);
PH = delayPhase(PH, freq, cdelay/1e3);
xplot(freq, PH, plotArgs, 'linestyle', 'none', 'markersize', 5);
xlabel('Frequency (kHz)')
ylabel('Phase (cycle)')
set(gca,'xtick',0:2:30);

function PG = local_trianplot(L1,L2,H1,H2, cdelay);
xplot([L1(1) L2(1)], [L1(2) L2(2)], 'k-');
xplot([H1(1) H2(1)], [H1(2) H2(2)], 'k-');
Lm = (L1+L2)/2; Hm = (H1+H2)/2;
DEL = @(P1,P2)-1e3*(P2(2)-P1(2))/(P2(1)-P1(1));
dL = cdelay+DEL(L1,L2)
dH = cdelay+DEL(H1,H2)
text(Lm(1), Lm(2)+0.07, sprintf('%d \\mus', round(dL)), ...
    'horizontalalign', 'right', 'fontweight', 'bold', 'backgroundcolor', 'w');
text(Hm(1), Hm(2)+0.1, sprintf('%d \\mus', round(dH)), ...
    'fontweight', 'bold', 'backgroundcolor', 'w');
PG = P2dB(dH/dL);
PG = 0.1*round(10*PG);
text(min(xlim), min(ylim)+0.067*diff(ylim),['+' num2str(cdelay) ' \mus']);

function local_markGain(ha, ab, Fx, PG);
axes(ha);
Gx = interp1(ab.Fprim/1e3, ab.Gain, Fx);
xplot([Fx Fx], Gx+0.5*PG*[-1 1], 'color', [1 0.8 0.8], 'linewidth', 3);

