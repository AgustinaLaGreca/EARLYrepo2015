function ThrCorrected = correct_ds(ds)
% CORRECT_DS is an auxiliary function that calculates the correction per
% frequency component for dataset ds for a THR curve.
% Introduced in order to fix the bug introduced
% in Apr 2018 (or earlier) in THR_curve.
%
% SpikeCrit can be obtained from an auxiliary file of the form THR_*.mat
% This variable is obtained in THRcurve for each THR, can be either custom
% manually introduced or calculated. It is used in the thrrec_procedure.
% It is not possible to use it to correct the curve but the code to
% retrieve is commented out.
% See THRcurve, THRrec_Geisler.
%


% Maximum time difference allowed between THR files saving times
er_time = 2/(24*60*60); % In days


experiment_dir = fileparts(filename(ds.ID.Experiment));
fname = [fileprefix(ds) '.EarlyDS'];

ThrDS = ds.Data.Thr;

% Unnecesary - Necessary to retrieve SpikeCrit stored in other file
% fdir = dir(fname);
% if ~isempty(fdir)
%     thr_timenum = fdir.datenum;
%     THR_files = dir(fullfile(experiment_dir, 'THR*'));
%     THR_times = [THR_files(:).datenum]';
%     idx = find(abs(THR_times - thr_timenum) < er_time, 1);
%     if ~isempty(idx)
%         dsTHR = load(fullfile(experiment_dir,THR_files(idx).name));
%     else; error('No THR_*.mat file found that matches this THR curve. SpikeCrit impossible to retrieve.');
%     end
%     % If the THR file found corresponds to the threshold dataset provided,
%     % retrieve SpikeCrit
%     ThrSpike = {dsTHR.S(:).Thr}';
%     
%     ThrSNan = ThrSpike;
%     ThrSNan(cellfun(@isempty,ThrSNan)) = {NaN};
%     idxThr = find(~cellfun(@isnan,ThrSNan));
%     Thr1 = cell2mat(ThrSpike(idxThr));
%     Thr2 = ThrDS(idxThr);
%     if isequal(Thr1,Thr2)
%         SpikeCrit = dsTHR.S(idxThr(1)).SpikeCrit;
%     else; error('No THR_*.mat file found that matches this THR curve. SpikeCrit impossible to retrieve.');
%     end
%     
% else; error(['No ' fname ' file found for the selected THR curve']);
% end

P = ds.Stim;

% Determine order
ifreqs = 1:numel(P.Freq);
if strcmpi(P.Order(1),'R')
    ifreqs = fliplr(ifreqs);
end

% Calculate # SPLs per Freq
NbSPL = size(P.LinAmp,1)/size(P.Freq,1);


[LinAmp_f, Attenuation_f] = local_getLinAmp(P,ifreqs, NbSPL, 'wrong');
[LinAmp_t, Attenuation_t] = local_getLinAmp(P,ifreqs, NbSPL, 'real');

% Calculate the correction between the LinAmp and Attenuation that were
% supposed to be used and the ones that were actually used
DC = zeros(size(ThrDS));
for n = 1:length(ThrDS)
    DC(n) = mean((20*log10(LinAmp_t{n})-Attenuation_t{n} - 20*log10(LinAmp_f{n})+Attenuation_f{n}),1);
end

ThrCorrected = ThrDS + DC;

end
    
function [LinAmp, Att, SPL] = local_getLinAmp(P, ifreqs, NbSPL, flag)

LinAmp = cell(1,numel(ifreqs));
Att = cell(1,numel(ifreqs));
SPL = cell(1,numel(ifreqs));

for ifreq=ifreqs
    
    ic = (1+(ifreq-1)*NbSPL:ifreq*NbSPL)';
    % Get correct LinAmp and Attenuations
    W               = P.Waveform;
    EXP             = P.Experiment;
    [MaxSPL, ~]     = maxSPL(W, EXP);
    SPLb            = P.SPLs;
    
    if strcmp(flag,'wrong')
        %  Bug introduced in April 2018 - see correct_ds
        SPL{ifreq}      = SPLb(SPLb < MaxSPL(ifreq));
        LinAmp{ifreq}   = P.LinAmp(1:length(SPL{ifreq}),:);
        Att{ifreq}      = P.Attenuations(1:length(SPL{ifreq}),:);
    else
        iSPL            = SPLb < MaxSPL(ifreq);  % eric 26/Mar/2019
        SPL{ifreq}      = SPLb(iSPL);            % eric 26/Mar/2019
        LinAmp{ifreq}   = P.LinAmp(ic(iSPL),:); % eric 26/Mar/2019
        Att{ifreq}      = P.Attenuations(ic(iSPL),:);   % eric 26/Mar/2019
    end
end
end