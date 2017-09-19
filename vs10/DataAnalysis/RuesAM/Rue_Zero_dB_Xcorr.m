function [S, allRho] = Rue_Zero_dB_Xcorr(FN1, FN2);
%  Rue_Zero_dB_Xcorr - correlation analysis of a recording pair
%    Single pair:
%    [S, allRho] = Rue_Zero_dB_Xcorr('p81213A', 'p81213B');
%    All true Pairs:
%     S = Rue_Zero_dB_Xcorr()
%    S is now a struct array.

if nargin<1, % all pairs
    RL = sort(RueList);
    RL = RL(1:2:end);
    for ipair=1:(numel(RL)/2),
        [FN1, FN2] = deal(RL{2*ipair-1}, RL{2*ipair})
        S(ipair) = Rue_Zero_dB_Xcorr(FN1, FN2);
    end
    return;
end;

% ====== single pair from here ========

% remove any trailing 'i' or 'c' from filenames
lastch = lower(FN1(end));
if isequal('i', lastch) || isequal('c', lastch),
    FN1 = FN1(1:end-1);
end
lastch = lower(FN2(end));
if isequal('i', lastch) || isequal('c', lastch),
    FN2 = FN2(1:end-1);
end
allRho = local_allRho(FN1,FN2);
idiag = 1:(size(allRho,1)+1):numel(allRho);
isecdiag = [0:(size(allRho,1)+1):numel(allRho) 2:(size(allRho,1)+1):numel(allRho)]; isecdiag(1)=[];
inondiag = setdiff(1:numel(allRho), idiag);
rho_d = allRho(idiag);
rho_2d = allRho(isecdiag);
rho_nd = allRho(inondiag);

mean_rho_diag = deciRound(mean(rho_d),3);
std_rho_diag = deciRound(std(rho_d),3);
mean_rho_nondiag = deciRound(mean(rho_nd),3);
std_rho_nondiag = deciRound(std(rho_nd),3);
mean_rho_secdiag = deciRound(mean(rho_2d),3);
std_rho_secdiag = deciRound(std(rho_2d),3);


ttest_diag001 = ttest(rho_d,0,0.01, 'both');
ttest_diag005 = ttest(rho_d,0,0.05, 'both');
ttest_nondiag001 = ttest(rho_nd,0,0.01, 'both');
ttest_nondiag005 = ttest(rho_nd,0,0.05, 'both');
ttest_secdiag001 = ttest(rho_2d,0,0.01, 'both');
ttest_secdiag005 = ttest(rho_2d,0,0.05, 'both');

S = CollectInStruct(FN1, FN2, '-', mean_rho_diag, std_rho_diag, ...
    mean_rho_nondiag, std_rho_nondiag, ...
    mean_rho_secdiag, std_rho_secdiag, ...
    '-', ...
    ttest_diag001, ttest_diag005, ...
    ttest_nondiag001, ttest_nondiag005, ...
    ttest_secdiag001, ttest_secdiag005);


function allRho = local_allRho(FN1,FN2);
PCache = {FN1,FN2};
allRho = getcache(mfilename, PCache);
if ~isempty(allRho),return; end
Nrep1 = RueNrep(FN1); Nrep2 = RueNrep(FN2);
allRho = nan(Nrep1,Nrep2);
for irep1=1:Nrep1,
    V1 = readRueRep(FN1, irep1);
    V1 = V1(:,505:end); % 0 dB only
    for irep2=1:Nrep2,
        V2 = readRueRep(FN2, irep2);
        V2 = V2(:,505:end); % 0 dB only
        allRho(irep1,irep2) = corr(V1(:), V2(:));
    end
end
putcache(mfilename, 100, PCache, allRho);


