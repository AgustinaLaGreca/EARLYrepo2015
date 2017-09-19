function S = tchebi(D, isweep, ichan, Nmax, Dwin, Vmin);
% tchebi - event extraction using tchebishev polynomials
% S = tchebi(D, isweep, ichan, Nmax, Dwin, Vmin);
%
%   See also getABFdata.

dt = D.dt_ms; % sample period in ms
Nwin = 1+round(Dwin/dt) % # samples in window
% generate chebi polynomials
t0 = linspace(-Dwin/2,Dwin/2,Nwin).';
W = sqrt(1-t0.^2)*sqrt(2/pi); % weighting function for chebis of 2nd kind
U = []; % unweighted chebis; U(:,n+1) is n-th-order chebi of 2nd kind
u = []; % weighted ones, i.e., u=W.*U
for iord =0:Nmax,
    qq = (iord+1)*hypergeom([-iord iord+2], 3/2, (1-t0)/2);
    U = [U, qq];
    u = [u, qq.*W];
end
%Normalize
NN = Nwin/2;
U = U/NN;
u = u/NN;
%plot(t0,u);

% convolve recording with chebis
x = D.AD(isweep,ichan).samples;
Y = [];
for iord=0:Nmax,
    Y = [Y, conv(x,u(end:-1:1,iord+1))];
end

% find "envelope"
%Env = sqrt(2*sum(Y.^2,2));
Env = 2*sqrt(conv(x.^2, W/sum(W)));
t = timeaxis(Env,dt)-Dwin/2;

% find "events"
[Tm, Em] = localmax(t,Env);
iok = find(Em>=Vmin);
Tm = Tm(iok);
Em = Em(iok);
Nevent = numel(Tm);

% reconstruct waveforms surrounding events using the chebis
Event = zeros(Nwin,Nevent);
Tcoef = zeros(Nevent, Nmax+1);
for ievent=1:Nevent,
    Tev = Tm(ievent); % event time in ms
    itime = 1+round((Tev+Dwin/2)/dt); % event-time index in Env & Y
    coef = Y(itime,:); % chebi coefficients of the event
    Event(:,ievent) = u*coef.'; % weighted sum of chebis
    Tcoef(ievent,:) = coef;
end
Event = NN*Event;

f1; set(gcf,'units', 'normalized', 'position', [0.047656 0.50391 0.4375 0.41016])
dplot(dt, x,'b')
xplot(t, Env,'r')
xplot(Tm,Em,'*k')
xlabel('time (ms)');
for ievent=1:Nevent,
    xplot(t0+Tm(ievent), Event(:,ievent), 'g');
end

FileName = D.Header.XX_name;
Param = CollectInStruct(FileName, isweep, ichan, Nmax, Dwin, Vmin);


DD=pdist(Tcoef); ZZ = linkage(DD);
iCluster = cluster(ZZ,'Cutoff',5,'Criterion','distance');
CLRs = get(gca,'colororder'); Nc = size(CLRs,1);

f2;  set(gcf,'units', 'normalized', 'position', [0.55938 0.49609 0.4375 0.41016])
for ii=1:max(iCluster), 
    CL = CLRs(1+rem(ii,Nc),:);
    ic = find(iCluster==ii);
    xplot(t0, Event(:,ic), 'color', CL);
    Cluster{ii} = ic;
end;

f3; set(gcf,'units', 'normalized', 'position', [0.36484 0.17125 0.3625 0.25125])
hist(DD,123);

disp('PCA!')
%[Coef, Score]=princomp(Tcoef(Cluster{end-1},:));



S = CollectInStruct(Param, t0,U,u,W,Event,Tcoef,iCluster, Cluster);
% plot(t,Y);
% xplot(t,sum(Y,2),'k', 'linewidth',2);
% xdplot(dt, x,'b:', 'linewidth',2)
% 
% 












