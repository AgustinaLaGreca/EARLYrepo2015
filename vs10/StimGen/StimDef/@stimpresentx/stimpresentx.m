function S = stimpresentx(Fsam, Nsample, Ncond, Nrep, iCond, iRep, X, Y, NsamBaseline, ReducedStorage);
% stimpresentx - constructor for stimpresentx objects
%    [Note: stimpresentx is updated version of stimpresent, including 
%     baseline specification . Added for backward compatibility.]
%    S = stimpresent(Fsam, Nsample, Ncond, Nrep, iCond, iRep, X, Y, BaselineDur)
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
% NsamBaseline: # samples of baseline intervals i.e., silent intervals
%             preceding and following the stimulus conditions, during which
%             recordings are made. Single values indicate identical pre-
%             and post-baseline sample counts. A 2-element array is
%             interpreted as [Pre, Post].
% ReducedStorage: criterion for reduced storage of presentations 
%                 to spare buffer space
%                 '' (unreduced) | 'nonzero' (only withold the nonzero
%                 segments of the first Waveform)
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
% ReducedStorage: criterion for reduced storage of pres
%   ------------------------------------
%         PreDur: duration (ms) of pre-stimulus baseline
%        PostDur: duration (ms) of post-stimulus baseline
%        NsamPre: # samples of pre-stimulus baseline
%       NsamPost: # samples of post-stimulus baseline
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
        iCond, iRep, NsamPres, SamOffset, PresOnset, PresDur, ReducedStorage, PreDur, PostDur, NsamPre, NsamPost, X, Y] = deal([]);
elseif isstruct(Fsam), % convert struct to stimpresentx
    S = Fsam;
    [Fsam, Npres, Ncond, Nrep, TotNsamPlay, TotDur, ...
        iCond, iRep, NsamPres, SamOffset, PresOnset, PresDur, PreDur, PostDur, NsamPre, NsamPost, X, Y] = deal([]);
    ReducedStorage = ''; % backwards compatible value
    Template = CollectInStruct(Fsam, Npres, Ncond, Nrep, TotNsamPlay, TotDur, '-', ...
    iCond, iRep, NsamPres, SamOffset, PresOnset, PresDur, ReducedStorage, '-', ...
    PreDur, PostDur,  NsamPre, NsamPost, '-', X, Y);
    S = union(Template, S);
    S = class(rmfield(S, 'dynamic'),mfilename,S.dynamic);
    return;
elseif isa(Fsam,mfilename) % stimpresentx(S): tautology
    S = struct(Fsam);
        [Fsam, Npres, Ncond, Nrep, TotNsamPlay, TotDur, ...
        iCond, iRep, NsamPres, SamOffset, PresOnset, PresDur, PreDur, PostDur, NsamPre, NsamPost, X, Y] = deal([]);
    ReducedStorage = ''; % backwards compatible value
    Template = CollectInStruct(Fsam, Npres, Ncond, Nrep, TotNsamPlay, TotDur, '-', ...
    iCond, iRep, NsamPres, SamOffset, PresOnset, PresDur, ReducedStorage, '-', ...
    PreDur, PostDur,  NsamPre, NsamPost, '-', X, Y);
    S = union(Template, S);
    S = class(rmfield(S, 'dynamic'),mfilename,S.dynamic);
    return;
else, % bookeeping
    % insert baseline intervals, if applicable
    [NsamPre, NsamPost] = DealElements(NsamBaseline([1 end]));
    [PreDur, PostDur] = DealElements([NsamPre, NsamPost]*1e3/Fsam);
    NstimCond = size(Nsample,1)-2; % # true stimulus conditions, not counting baselines
    if NsamPre>0,
        iCond = [NstimCond+1; iCond]; % condition # Ncond+1 is pre-stim baseline by def
        iRep = [1; iRep]; 
    end
    if NsamPost>0,
        iCond = [iCond; NstimCond+2]; % condition # Ncond+2 is post-stim baseline by def
        iRep = [iRep; 1]; 
    end
    if ~isequal(Nsample(:,1), Nsample(:,end)),
        error('Left & Right-channel waveforms contain unequal number of samples.');
    end
    Nsample = Nsample(:,1);
    Npres = numel(iCond); % number of individual presentations, counting reps

    NsamPres = Nsample(iCond); % sample count per presentation
    TotNsamPlay = sum(NsamPres);
    TotDur = 1e3*TotNsamPlay/Fsam; % total stimulus duration in ms
    % Construct a fcn that converts sample indices to presentation counts
    SamOffset = [0; cumsum(NsamPres)]; % sample indices corresponding to presentation boundaries
    PresOnset = 1e3*SamOffset/Fsam; % onsets [ms] of presentations
    PresDur = diff(PresOnset); % exact durations of the presentations
    PresDur = deciRound(PresDur, 10); % round to 10 decimals to equalize small rounding errors
    
end

% store all info in Presentation field of P
S = CollectInStruct(Fsam, Npres, Ncond, Nrep, TotNsamPlay, TotDur, '-', ...
    iCond, iRep, NsamPres, SamOffset, PresOnset, PresDur, ReducedStorage, '-', ...
    PreDur, PostDur,  NsamPre, NsamPost, '-', X, Y);
S = class(S,mfilename,dynamic);













