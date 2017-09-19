function S = RueXXXcorr(FN1, FN2);
% RueXXXcorr - the ultimate xcor for Rue's audiogram data
%   S = RueXXXcorr(FN1, FN2) returns struct array S with fields
%        FN1: cell #1
%        FN2: cell #2
%     isPair: true for a cell pair
%     isSame: true when FN1 & FN2 are the same
%
%       freq: tone frequency in kHz
%        SPL: sound intensity in dB SPL
%      Nrep1: # repetitions for cell 1
%      Nrep2: # repetitions for cell 2
%
%      mean1: mean values of waveforms of reps # 1..Nrep1 of cell 1
%      mean2: mean values of waveforms of reps # 1..Nrep2 of cell 2
%      var1: variances of waveforms of reps # 1..Nrep1 of cell 1
%      var2: variances of waveforms of reps # 1..Nrep2 of cell 2
%
%       rho: normalized Xcorr values in Nrep1xNrep2 matrix
%      urho: unnormalized  Xcorr values in Nrep1xNrep2 matrix
%    meanrho: mean across all rep-pairs of rho
%    meanrho_nd: mean across all nondiagonal rep-pairs of rho
%
%     MWvar1: variance of waveform average across reps of cell 1
%     MWvar2: variance of waveform average across reps of cell 2
%      MWrho: normalized Xcorr between mean waveforms
%     MWurho: unnormalized Xcorr between mean waveforms

if nargin<2,
    FN2 = FN1;
elseif isequal('B', FN2),
    FN2 = strrep(FN1, 'A', 'B');
end

% sort F1 and F2 
FN = sort({FN1 FN2});
FN1 = FN{1}; FN2 = FN{2};

% try to retrieve S from cache
CFN = mfilename;
CPar = {FN1, FN2, 1};
S = getcache(CFN, CPar);
if ~isempty(S), return; end

disp(['=== ' FN1, '  ' FN2 '===']);
% compute all fields of S
isPair = ~isequal(FN1, FN2) && (...
    isequal(strrep(FN1,'A', 'B'), FN2) ...
    || isequal(strrep(FN1,'B', 'A'), FN2)   );
isSame = isequal(FN1, FN2);

allSPLs = 80:-10:0; % dB SPL
allFreqs = 2.^((0:28)/5); % kHz

S = [];
for iSPL=1:numel(allSPLs),
    SPL = allSPLs(iSPL);
    disp(['  ' num2str(SPL) ' dB'])
    for ifreq=1:numel(allFreqs),
        freq = allFreqs(ifreq);
        [Y1, Nrep1, dt] = readRueCond(FN1, iSPL, ifreq);
        [Y2, Nrep2, dt] = readRueCond(FN2, iSPL, ifreq);
        % subplot(2,1,1); dplot(dt, Y1); xdplot(dt, mean(Y1,2),'k', 'linewidth', 4);
        % subplot(2,1,2); dplot(dt, Y2); xdplot(dt, mean(Y2,2),'k',
        % 'linewidth', 4);
        mean1 = mean(Y1); % mean of each rep
        mean2 = mean(Y2); % mean of each rep
        var1 = var(Y1); % mean of each rep
        var2 = var(Y2); % mean of each rep
        rho = corr(Y1,Y2); % all crosscorrs across columns of Y1 and Y2
        [v1, v2] = SameSize(var1(:), var2(:).'); 
        urho = sqrt(v1.*v2).*rho;
        meanrho = mean(rho(:));
        rr = rho; for ii=1:min(Nrep1,Nrep2), rr(ii,ii)=nan; end
        meanrho_nd = mean(denan(rr(:)));
        % metrics MW.. based on mean waveforms
        ym1 = mean(Y1,2); ym2 = mean(Y2,2);
        MWvar1 = var(ym1);
        MWvar2 = var(ym2);
        MWrho = corr(ym1, ym2); % corr of means across reps
        MWurho  = std(ym1)*std(ym2)*MWrho; % unnormalized version of rhom
        S = [S, ...
            CollectInStruct(...
            FN1, FN2, isPair, isSame, '-', ...
            freq, SPL, Nrep1, Nrep2, '-', ...
            mean1, mean2, var1, var2, rho, urho, meanrho, meanrho_nd, '-', ...
            MWvar1, MWvar2, MWrho, MWurho);
            ];
    end
end
putcache(CFN, 1e4, CPar, S);


% ================================







