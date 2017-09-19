function S = stimpresent(Fsam, Nsample, Ncond, Nrep, iCond, iRep, X, Y);
% stimpresent - constructor for stimpresent objects
%    S = stimpresent(Fsam, Nsample, Ncond, Nrep, iCond, iRep, X, Y)
%    returns a StimPresent object containing  the details of a stimulus
%    presentation. Inputs:
%      Fsam: exact sample rate in Hz
%   Nsample: Nsample(j,ichan) is # samples played over DA channel ichan
%            for stimulus condition (j).
%     Ncond: number of different stimulus conditions sorted in "forward
%            direction" according to the stimGUI
%      Nrep: number of repetitions per stim conditions
%     iCond: iCond(k) is the index of the condition played during the
%            k-th presentation
%      iRep: iRep(k) is the repetition count of the k-th presentation
%      X, Y: structs containing details on varied parameters. If only one
%            stimulus parameter is varied, it is described in X.
%
%   Note: a "condition" is a stimulus corresponding to a certain value of
%   the stimulus parameter(s), e.g., 1500 Hz carrier frequency. A
%   "presentation" is an individual repetition of a condition. Therefore,
%   Npres == Ncond*Nrep.
%
%   S is a dynamic object (see Dynamic).
%   struct(S) has the containing the fields (see also S.help)
%           Fsam:  sample rate [Hz]
%          Npres: number of stimulus presentations, reps counted
%          Ncond: number of different conditions
%           Nrep: number of repetitions per condition
%    TotNsamPlay: total # samples in stimulus
%         TotDur: total stim duration [ms]
%    -----------------------------------
%          iCond: index array; iCond(k)==j -> k-th pres plays j-th condition
%           iRep: index array; iRep(k)==j -> k-th pres is the j-th rep of that condition
%       NsamPres: sample counts; NsamPres(k)=n -> k-th pres plays n samples
%      SamOffset: sample offsets of pres since start of stim
%      PresOnset: time offsets [ms] of pres since start of stim
%        PresDur: pres durations [ms]
%    -----------------------------------
%              X: details on 1st varied stimulus param [struct]
%              Y: details on 2nd varied stimulus param [struct]
%
%    StimPresent called with no input arguments returns a void stimpresent
%    object.
%
%    StimPresent is typically called by SortConditions.
%
%    See also SortConditions, action, dynamic.

if nargin<1, % void object
    [Fsam, Npres, Ncond, Nrep, TotNsamPlay, TotDur, ...
        iCond, iRep, NsamPres, SamOffset, PresOnset, PresDur, X, Y] = deal([]);
else,
    if ~isequal(Nsample(:,1), Nsample(:,end)),
        error('Left & Right-channel waveforms contain unequal number of samples.');
    end
    Nsample = Nsample(:,1);
    Npres = numel(iCond);; % number of individual presentations, counting reps

    NsamPres = Nsample(iCond); % sample count per presentation
    TotNsamPlay = sum(NsamPres);
    TotDur = 1e3*TotNsamPlay/Fsam; % total stimulus duration in ms
    % Construct a fcn that converts sample indices to presentation counts
    SamOffset = [0; cumsum(NsamPres)]; % sample indices corresponding to presentation boundaries
    PresOnset = 1e3*SamOffset/Fsam; % onsets [ms] of presentations
    PresDur = diff(PresOnset); % exact durations of the presentations
    RemainingDur = TotDur-PresOnset;
end

% store all info in Presentation field of P
S = CollectInStruct(Fsam, Npres, Ncond, Nrep, TotNsamPlay, TotDur, '-', ...
    iCond, iRep, NsamPres, SamOffset, PresOnset, PresDur, '-', X, Y);
S = class(S,mfilename,dynamic);













