function [SPT, Nspike] = spiketimes(DS, Chan, flag, minISI);
% Dataset/SpikeTimes - sorted spike times of dataset
%    SPT = SpikeTimes(DS) returns the spike times of dataset DS in a cell 
%    array SPT, with SPT{iCond, iRep} an array holding the the spike 
%    arrival times [ms] re the onset of the stimulus presentation in 
%    response to repetition iRep of stimulus condition iCond. 
%
%    NOTE: If offline spike times exist for dataset DS, these will overrule
%    any spikes saved during the recording. See dataset/HasOfflineSpiketimes.
%
%    [SPT, Nspike] = SpikeTimes(DS) also returns an numerical array Nspike
%    with Nspike(iCond,iRep)==numel(SPT{iCond,iRep}), that is, the spike
%    count in response to stimulus presentation (iCond,iRep).
%
%    SpikeTimes(D, Chan) extracts the spiketimes from the digital channel
%    specified in Chan. See Digichan for naming conventions.
%    The default Chan is 1 or, equivalently,'RX6_digital_1'.
%
%    SpikeTimes(D, Chan, 'no-unwarp') cancels time-unwarping of the spike 
%    times according to the Warp value of each stimulus condition. The 
%    default is that warped conditions are automatically unwarped.
%    SpikeTimes(D, Chan, '') does apply unwarping and is the same 
%    as SpikeTimes(D, Chan). See stimDefMask. 
%
%    SpikeTimes(D, Chan, flag, minISI) removes any spikes following their
%    predecessors within minISI ms. Default minISI = 0 ms, i.e., no
%    elimination.
%
%    See also Dataset/SPTbaseline, Dataset/dotraster, Dataset/AnWin.        
[Chan, flag, minISI] = arginDefaults('Chan/flag/minISI',1, '', 0);

if isvoid(DS),
    error('Void dataset.');
end
Unwarp = true;
if isequal('', flag),
elseif isequal('no-unwarp', flag),
    Unwarp = false;
else,
    error('invalid flag. Must be empty or ''no-unwarp''.');
end

PRES = DS.Stim.Presentation;

if HasOfflineSpikes(DS), % use previously saved offline spike times
    SPT = getfield(load(OfflineSpikeFilename(DS), '-mat'), 'SPT');
else, % use events stored online contained in DS
    % for each event, retrieve most recent stimulus onset & its index
    ET = digichan(DS, Chan);
    spt = eventtimes(ET); % spike times re grand stimulus onset
    [Onsets, ipres] = NthFloor(spt, PRES.PresOnset(2:end)); % per spike, the most recent onset (ignore 1 = baseline)
    spt = spt-Onsets; % spike times re "their own" stim onsets
    % put spike times in NcondxNrep cell matrix as described in help text
    SPT = cell(PRES.Ncond, PRES.Nrep);
    
    % Check if the stutter was on => Ignore the first presentation after
    % the baseline; not needed for presentation
    if isfield(DS.Stim,'Stutter')
        StutterOn = 1;
    else
        StutterOn = 0;
    end

    %%%%%%%%%%%%%%%% Code Adriaan to avoid continious pointer dereferencing in for loop 
    icond_arr = DS.Stim.Presentation.iCond;
    irep_arr = DS.Stim.Presentation.iRep;

    for ip=1:PRES.Ncond*PRES.Nrep,       
        if ~StutterOn
            icond = icond_arr(ip+1); % condition index of presentation ip; +1 to account for baseline
            irep = irep_arr(ip+1); % repetition count of presentation ip; +1 to account for baseline
            SPT{icond,irep} = spt(ipres==ip);  
        else
            icond = icond_arr(ip+2); % condition index of presentation ip; +2 to account for baseline and stutter
            irep = irep_arr(ip+1); % repetition count of presentation ip; +1 to account for baseline
            SPT{icond,irep} = spt(ipres==(ip+1));
        end
    end

    clear icond_arr irep_arr;
end
    
Nspike = cellfun('length',SPT);

if Unwarp,
    St = DS.Stim;
    if isfield(St, 'Warp'), 
        for icond = 1:PRES.Ncond,
            SPT(icond,:) = local_dewarp(SPT, icond, St.Warp);
        end
    end
end

if minISI>0, % eliminate spike "echos" following previous spike too fast
    for ii=1:numel(SPT),
        ifast = 1+find(diff(SPT{ii})<minISI);
        SPT{ii}(ifast) = [];
    end
end

%======================================

function  spt = local_dewarp(SPT, icond, Warp);
spt =  SPT(icond,:);
warpfac = 2^Warp(icond);
warpfac = SameSize({warpfac}, spt);
spt = cellfun(@times, spt, warpfac, 'UniformOutput', false);






