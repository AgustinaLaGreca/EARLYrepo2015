function [S, Details]=deristats(D, tau, ClickDownSlopeAP, ClickDownSlopeEPSP, ClickUpSlopeAP);
% deristats - view statistics of time-derivative of recordings 
%   deristats(D, tau) 
%      D: ABF file in TK format (see readTKABF)
%    tau: width (ms) of smooting hann window. Default 0.2 ms.

if nargin<3, ClickDownSlopeAP=[]; end % default: ask user to click in graph
if nargin<4, ClickDownSlopeEPSP=[]; end % default: ask user to click in graph
if nargin<5, ClickUpSlopeAP=[]; end % default: ask user to click in graph

MaxIntPreAP = 0.8; % ms max interval preceding AP peak in which to look for event onset

Npart = 5;
CM = colormap; Ncol = size(CM,1);
ico = round(linspace(1,Ncol,Npart));
LiCo = CM(ico,:);
MinFallDur = 0.2; % ms min duration of fall portions to be analyzed
MinVfall = 0.3; % mV min drop during falling portions

if nargin<2, tau = 0.2; end % ms default windiw width

sam = catsam(D); % whole recording
%sam = catsam(D,2); '=======DEBUG========'
MeanSam = mean(sam);
MedSam = median(sam); % median value of samples, serving as "baseline"
IqrSam = iqr(sam); % inter quartile range, serving as "natural spread"
%sam = cat(1,D.AD(2,1).samples); % debug

Nsam = numel(sam);

% smoothing
Nwin = 1+2*round((tau/D.dt_ms-1)/2); % # samples of hann window
hWin = hann(Nwin); hWin = hWin/sum(hWin) % normalized hann window
smsam = conv(sam, hWin); % smoothed version of the samples
Nhalf = round((Nwin-1)/2); fakeWin = 0*hWin + [zeros(Nhalf,1); 1; zeros(Nhalf,1)];
sam = conv(sam, fakeWin); % unsmoothed version of the samples having same size as smsam

dsam = diff(smsam)/D.dt_ms; % time derivative V/S
[dmin, dmax] = minmax(dsam);
Edges = linspace(dmin, dmax*1.0001, 300)';
N = hist(dsam, Edges); N = N(1:end-1); N = N(:);
BinCenters = (Edges(1:end-1)+Edges(2:end))/2;
xplot(BinCenters, N);
T = timeaxis(sam,D.dt_ms)-Nhalf*D.dt_ms;


% if isnan(Ymin), Ymin=0.05*max(N); end
% Yfit = N(N>Ymin);
% Bfit = BinCenters(N>Ymin);
% FM = fit(Bfit, Yfit, 'gauss1');
% Ym = feval(FM, Bfit);
% 
% xplot(Bfit, Ym, 'r');
% xplot(Bfit, Yfit-Ym, 'k:');

E2 = linspace(dmin, dmax, Npart+1);
for ipart=1:Npart,
    ihit = find(betwixt(dsam, E2(ipart+[0 1])));
    Ip{ipart} = ihit(betwixt(ihit, [0 Nsam+1]));
end


% analyze time constants of rising and falling portions
istartrise = 1+find((dsam(1:end-1)<=0) & (dsam(2:end)>0)); 
istartfall = 1+find((dsam(1:end-1)>0) & (dsam(2:end)<=0)); 
% fenceplot(T(istartrise), ylim, 'g');
% fenceplot(T(istartfall), ylim, 'r');

% find falling portions
iturn = sort([istartrise(:); istartfall(:)]);
ifall = iturn; 
if ifall(1)==istartrise(1), ifall(1)=[]; end % make sure first turn is downward
if rem(numel(ifall),2)>0, ifall(end)=[]; end % make sure last turn is upward
ifall = reshape(ifall, [2,round(numel(ifall)/2)]);
% ditto rising portions
irise = iturn;
if irise(1)==istartfall(1), irise(1)=[]; end % make sure first turn is upward
if rem(numel(irise),2)>0, irise(end)=[]; end % make sure last turn is downward
irise = reshape(irise, [2,round(numel(irise)/2)]);
FirstIsRising = irise(1)<ifall(1);
% fenceplot(T(iturn(1,:)), ylim, 'g');
% fenceplot(T(iturn(2,:)), ylim, 'r');
;


%==== use user input for deparating APs, EPSPs, noise
figure; set(gcf,'units', 'normalized', 'position', [0.357 0.14 0.361 0.76]);
% APs
samFall = sam(ifall); 
MagnFall = -diff(samFall,1,1); % size of MagnFall (mV) during fall portions
Tfall = T(ifall); % time (ms) at turning points
FallDur = diff(Tfall,1,1); % durations of falling portions
plot(FallDur, MagnFall, 'r.');
if isempty(ClickDownSlopeAP), % prompt user
    uiwait(msgbox('Left-click twice on graph to separate APs from the rest.', 'user input'));
    ClickDownSlopeAP=ginput(2);
end
X=ClickDownSlopeAP(:,1); Y=ClickDownSlopeAP(:,2);
xplot(X,Y,'k-');
Normal = [Y(2)-Y(1), X(1)-X(2)];
Normal = -Normal*sign(Normal(2));
Const = Normal*[X(:)'; Y(:)'];
isAP = (Normal*[FallDur; MagnFall])<Const(1);
xplot(FallDur(isAP), MagnFall(isAP), 'b.');
% EPSPs
dsize(MagnFall, isAP, MagnFall(isAP));
ylim([0 1.1*max(max(MagnFall(~isAP)))]);
if isempty(ClickDownSlopeEPSP), % prompt user
    uiwait(msgbox('Left-click twice on graph to separate EPSPs from noise.', 'user input'));
    ClickDownSlopeEPSP=ginput(2);
end
X=ClickDownSlopeEPSP(:,1); Y=ClickDownSlopeEPSP(:,2);
xplot(X,Y,'k-');
Normal = [Y(2)-Y(1), X(1)-X(2)];
Normal = -Normal*sign(Normal(2));
Const = Normal*[X(:)'; Y(:)'];
isEvent = (Normal*[FallDur; MagnFall])<Const(1);
xplot(FallDur(~isEvent), MagnFall(~isEvent), '.', 'color', 0.8*[1 1 1]);
isEPSP = isEvent & (~isAP);
ylim auto;

%==== use user input for separating APs, EPSPs, noise
figure; set(gcf,'units', 'normalized', 'position', [0.357 0.14 0.361 0.76]);
% APs
MagnRise = diff(sam(irise),1,1); % size of rise (mV) during rise portions
Trise = T(irise); % time (ms) at turning points
RiseDur = diff(Trise,1,1); % durations of riseing portions
MeanUpSlope = MagnRise./RiseDur;
plot(RiseDur, MagnRise, 'r.');
if isempty(ClickUpSlopeAP), % prompt user
    uiwait(msgbox('Left-click twice on graph to separate steepest stuff from the rest.', 'user input'));
    ClickUpSlopeAP=ginput(2);
end
X=ClickUpSlopeAP(:,1); Y=ClickUpSlopeAP(:,2);
xplot(X,Y,'k-');
Normal = [Y(2)-Y(1), X(1)-X(2)];
Normal = -Normal*sign(Normal(2));
Const = Normal*[X(:)'; Y(:)'];
isSteep = (Normal*[RiseDur; MagnRise])<Const(1);
xplot(RiseDur(isSteep), MagnRise(isSteep), 'b.');
% EPSPs
% dsize(MagnRise, isAP, MagnRise(isAP));
% ylim([0 1.1*max(max(MagnRise(~isAP)))]);
% uiwait(msgbox('Left-click twice on graph to separate EPSPs from noise.', 'user input'));
% qq=ginput(2); X=qq(:,1); Y=qq(:,2);
% xplot(X,Y,'k-');
% Normal = [Y(2)-Y(1), X(1)-X(2)];
% Normal = -Normal*sign(Normal(2));
% Const = Normal*[X(:)'; Y(:)'];
% isEvent = (Normal*[RiseDur; MagnRise])<Const(1);
% xplot(RiseDur(~isEvent), MagnRise(~isEvent), '.', 'color', 0.8*[1 1 1]);
% isEPSP = isEvent & (~isAP);
ylim auto;


% plot recording trace using color code for derivative values
hft = figure; set(gcf,'units', 'normalized', 'position', [0.0203 0.346 0.976 0.525])
plot(T,sam, 'k-');
xplot(T,smsam, '-', 'color', 0.7*[1 1 1]);
for ipart=1:Npart,
    ip = Ip{ipart};
    xplot(T(ip), sam(ip), 'color', LiCo(ipart,:), 'linestyle', 'none', 'marker', '.');
end
ylim(ylim);
xplot(xlim, MedSam+[0 0],'k:')
Ybar = 0; %mean(ylim)+0.45*diff(ylim);
dsize(Tfall(:, isAP), Ybar+0*Tfall(:, isAP) )
xplot(Tfall(:, isAP), Ybar+0*Tfall(:, isAP), 'r', 'linewidth', 4);
xplot(Tfall(:, isEPSP), Ybar+0*Tfall(:, isEPSP), 'b', 'linewidth', 3);

xplot(Trise(:, isSteep), Ybar+0*Tfall(:, isSteep), 'm', 'linewidth', 4);
%xplot(Tfall(:, isEPSP), Ybar+0*Tfall(:, isEPSP), 'b', 'linewidth', 3);

ylim(ylim);
DS = sgsrdataset(D.ExpID, D.RecID);
IDstr = [D.ExpID '/' D.RecID '/' num2str(D.icond) ' ---- ' num2str(DS.Xval(D.icond),3) ' ' DS.Xunit];
title(IDstr);

% analyze intervals just preceding AP down slopes
NsamPre = round(MaxIntPreAP/D.dt_ms);
irangePre = -NsamPre:-1;
iAP = find(isAP); 
ThrOnset = MedSam + 1.5*IqrSam; % threshold for baseline-like value
ionsetAP = nan(1,numel(iAP));
isteepEPSP_AP = nan(1,numel(iAP));
for ia=1:numel(iAP),
    isamPeak = ifall(1,iAP(ia)); % sample idx @ AP peak 
    isampre = isamPeak+irangePre; isampre(isampre<1)=[];
    presam  = smsam(isampre); % smoothed samples preceding  AP peak
    if ~isempty(presam), 
        ibase = isampre(find(presam<ThrOnset, 1, 'last')); % latest baseline-like sam before AP
        if isempty(ibase), continue; end
        io = ifall(2, find(ifall(2,:)<=ibase, 1, 'last')); % end of most recent fall preceding this baseline sample
        if isempty(io), continue; end
        if abs(sam(io)-MedSam)>2*IqrSam, % this rising phase does not start from ~ baseline;  ...
            ionsetAP(ia) = ibase; % ... take ibase as starting point
        else,
            ionsetAP(ia) = io;
        end
        isamrise = io:isamPeak-1;
        % determine steepest part of EPSP
        isteep = localmax(isamrise, dsam(isamrise));
        if isempty(isteep), continue; end
        if numel(isteep)>1,
            isteepEPSP_AP(ia) = isteep(end-1);
        else,
            isteepEPSP_AP(ia) = isteep(1);
        end
    end
end
xplot(T(denan(ionsetAP)), sam(denan(ionsetAP)), 'd', 'color', 0.3*[1 1 1]);
%xplot(T(isteepEPSP_AP), sam(isteepEPSP_AP), '^', 'color', [1 0 0]);

% derivative plot
% figure;
% plot(T(1:end-1), dsam);
% xplot(T(ionsetAP), sam(ionsetAP), 'd', 'color', 0.3*[1 1 1]);
% xplot(Tfall(:, isAP), Ybar+0*Tfall(:, isAP), 'r', 'linewidth', 4);
% xplot(Tfall(:, isEPSP), Ybar+0*Tfall(:, isEPSP), 'b', 'linewidth', 3);
% xplot(Trise(:, isSteep), Ybar+0*Tfall(:, isSteep), 'm', 'linewidth', 4);

%MinFallDur, MinVfall
ExpID = D.ExpID; RecID = D.RecID; icond = D.icond; 
S = CollectInStruct(ExpID ,RecID, icond, '-', ...
    tau, ClickDownSlopeAP, ClickDownSlopeEPSP, ClickUpSlopeAP, '-', ...
    MeanSam, MedSam, IqrSam, '-', ...
    Tfall, FallDur, MagnFall, samFall, isAP, isEPSP, FirstIsRising, '-', ...
    ionsetAP, isteepEPSP_AP);

Details = CollectInStruct(sam, dsam, smsam, T);
% g0=fittype('a*exp(-(x/c).^2)');
% Y=fit(BinCenters(:), N(:), g0, 'MaxFunEvals', 1e5);
% Nfit = feval(Y,BinCenters);
% xplot(BinCenters, Nfit, 'r');
% xplot(BinCenters, N-Nfit, 'k:')