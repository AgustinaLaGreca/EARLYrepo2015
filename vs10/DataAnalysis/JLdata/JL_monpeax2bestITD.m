function M= JL_monpeax2bestITD(Uidx);
%  JL_monpeax2bestITD - prdecit best ITD from mon peaks of bin recording
%      S = JL_monpeax2bestITD(Uidx);
%      Database in D:\processed_data\JL\monpeax2bestITD\monpeax2bestITD.mat
%      To build it:
%        qq = JL_rec_select;
%        qq = qq([qq.CanBeUsed]);
%        for ii=1:numel(qq); 
%           JL_monpeax2bestITD(qq(ii).UniqueRecordingIndex);
%           disp(ii);
%        end
%      To retrieve it:
%      DBFN = fullfile('D:\processed_data\JL\monpeax2bestITD', 'JL_monpeax2bestITD.dbase');
%      qqq=retrieve_dbase(DBFN)

S = JLdatastruct(Uidx); % specific info on stimulus etc
CFN = fullfile('D:\processed_data\JL\monpeax2bestITD', [mfilename '_' num2str(S.UniqueSeriesIndex)]);
[M, CFN, CP] = getcache(CFN, Uidx);

if isempty(M), % compute it and cache it
    D = JLdbase(Uidx); % general info on the recording & stimulus
    W = JLwaveforms(Uidx); % waveforms and duration params
    T = JLcycleStats(Uidx); % cycle-based means & variances
    SP = JLspikeStats(Uidx);
    %[T, Tstats] = JLanova2(Uidx);

    Ph_ipsiPk = local_phipeak(W.IpsiMeanrec);
    Ph_contraPk = local_phipeak(W.ContraMeanrec);
    predBestIPD = mod(Ph_contraPk-Ph_ipsiPk+0.5,1)-0.5;

    Nspikes = sum(W.APinStim);
    measBestIPD = mod(cangle(SP.VS_B)+0.5,1)-0.5;
    VS_BestIPD = abs(SP.VS_B);
    alphaBestIPD = SP.vs_alpha_B;
    Freq = S.Freq1;
    SPL = max([S.SPL1 S.SPL2]);
    StimType = S.StimType;
    UniqueRecordingIndex = Uidx;
    M = CollectInStruct(UniqueRecordingIndex, '-monpeax2bestITD',  ...
        SPL, Freq, StimType, Ph_ipsiPk, Ph_contraPk, predBestIPD, ...
        '-binSpikeTiming', Nspikes, measBestIPD, VS_BestIPD, alphaBestIPD);
    putcache(CFN, 50, Uidx, M);
    % add to database
    DBFN = fullfile('D:\processed_data\JL\monpeax2bestITD', [mfilename '.dbase']);
    init_dbase(DBFN, 'UniqueRecordingIndex', 'onlyifnew');
    Add2dbase(DBFN, M);
end



% ====================================
function ph = local_phipeak(Wv);
% peak of one cycle
dphi = 1/numel(Wv);
[ph, dum, isort] = localmax(dphi, Wv);
if isempty(isort),
    ph = nan;
else,
    ph = ph(isort(1));
end








