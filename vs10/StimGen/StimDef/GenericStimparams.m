function P = GenericStimparams(S, ParamName, iCond, ichan);
%  GenericStimparams - generic set of stimulus parameters
%     GenericStimparams(S), where S is the stimulus-defining struct of a
%     stimulus protocol (as stored in the Stim field of dataset objects),
%     returns a struct whose fields are a fixed, generic set of stimulus 
%     parameters. This fixed set of parameters has to be supplied for each
%     stimulus protocol FOO by the function stimparamsFOO.
%     All parameters that may vary across conditions contain Ncond 
%     rows. All parameters that may differ across D/A channels must contain
%     two columns. An error is thrown if the information supplied by
%     stimparamsFOO is not compatible with these specs, or lacks any of the
%     parameters listed in the table below.
%
%     GenericStimparams(S, 'ABC') only returns parameter ABC.
%
%     GenericStimparams(S, 'ABC', iCond) only returns the value of ABC for
%     condition(s) # iCond, that is, ABC(iCond,:). 
%
%     GenericStimparams(S, 'ABC', iCond, ichan) only returns the value of 
%     ABC for condition(s) # iCond and channel(s) ichan, that is, ABC(iCond, ichan). 
%     Specify iCond==0 if you want all conditions, but a only single channel.
%
%     GenericStimparams(D, ...), where D is a dataset object, is equivalent
%     to GenericStimparams(D.Stim, ...).
%
%     The parameters are described in the following table. ["pres" means
%     single stimulus presentation (=repetition); "rec" means recorded 
%     waveform for a single pres; the 3D arrays Fcar and CarStartPhase must
%     be squeezed into 2 dimensions if Ntone==1].
%     
%               Name    Type     Unit     Size        Meaning
%     ------------------------------------------------------
%     ------------------------------------------------------
%            StimType   Char      -         -      stimulus type (protocol)
%               Ncond   Double    -       1 x 1    # stim conditions
%               Nrep    Double   -        1 x 1    # repetitions per condition 
%               Ntone   Double   -        1 x 1    # simultaneous tones in one pres 
%     ------------------------------------------------------
%      PreBaselineDur   double   ms        1 x 1   silence preceding first pres
%     PostBaselineDur   double   ms        1 x 1   silence following last pres  
%                 ISI   double   ms    Ncond x 1   interval btw pres onset and next pres onset 
%            BurstDur   double   ms    Ncond x 2   duration of sound 
%          OnsetDelay   double   ms    Ncond x 1   delay from rec start to pres onset (excluding ITDs)
%             RiseDur   double   ms    Ncond x 2   duration of rising phase 
%             FallDur   double   ms    Ncond x 2   duration of falling phase 
%                 ITD   double   ms    Ncond x 1   interaural time delay
%             ITDtype:  char      -          -     type of ITD (e.g. 'ongoing')
%      TimeWarpFactor   double   ms    Ncond x 1   factor of time scaling (usually 1, but see MASK)
%     ------------------------------------------------------
%                Fsam   double   Hz       1 x 1      D/A sample rate
%                Fcar   double   Hz  Ncond x 2 x Ntone  carrier frequency 
%                Fmod   double   Hz     Ncond x 2    modulation frequency
%           LowCutoff   double   Hz     Ncond x 2    low cutoff frequency
%          HighCutoff   double   Hz     Ncond x 2    high cutoff frequency
%      FreqWarpFactor   double   Hz     Ncond x 1    reciprocal of TimeWarpFactor  
%     ------------------------------------------------------
%       CarStartPhase   double cycle Ncond x 2 x Ntone  carrier starting phase
%       ModStartPhase   double   cycle  Ncond x 2    modulation starting phase 
%            ModTheta   double   cycle  Ncond x 2    modulation angle (AM,QFM,mixed)
%            ModDepth   double   %      Ncond x 2    modulation depth
%     ------------------------------------------------------
%                 SPL   double  dB SPL  Ncond x 2    sound level
%             SPLtype   char      -        -         type of SPL, e.g. 'per tone'  
%                 DAC:  char      -        -         DA channels used (Both|Left|Right)
%     ------------------------------------------------------
%
%    NOTE: OnsetDelay must correspond to an integer number of samples, as
%    it is realized by prepending zero samples to the waveforms.
%
%    See also stimparamsFS.

[ParamName, iCond, ichan] = arginDefaults('ParamName, iCond, ichan', [],[],[]);

if isa(S, 'dataset'), S = S.Stim; end % see help text
% either extract info already contained in S, or delegate to stimparamsXXX
P = [];
if isfield(S,'GenericStimParams'), % already evaluated - check validity
    P = S.GenericStimParams;
    P = local_test(P); % if P fails the test, local_test returns [] 
end

if isempty(P), 
    spf = ['stimparams' S.StimType]; % name of stimulus-specific stimparams file ("stimParamsFOO")
    if isfield(S, 'GenericParamsCall'), % first choice: get from feval-ing P.GenericParamsCall.
        gpc = S.GenericParamsCall;
        gpc{2} = S;
        P = feval(gpc{:});
    elseif exist(spf, 'file'), % Second choice: stimParamsFOO file
        P = feval(spf, S);
        P.CreatedBy = spf;
    else,
        warning(['Cannot find parameter extractor ' spf '. Provide it!']);
        return;
    end

    [P, Mess] = local_test(P); 
    error(Mess);
end
%  further selection only if requested
if isempty(ParamName), return; end
P = P.(ParamName);
if isempty(iCond), return; end
if ~isequal(0, iCond),
    P = P(iCond,:,:); % 3rd dim may apply to Fcar, but doesn't do any harm for others
end
if isempty(ichan), return; end
P = P(:,ichan,:); % 3rd dim may apply to Fcar, but doesn't do any harm for others


%================================
function [P, Mess]=local_test(p);
Mess = ''; 
P = [];
% special case
if isequal('waiver',p), return; end;
% test for required field names
RFNS = {'StimType'    'Ncond'     'PreBaselineDur'    'PostBaselineDur'    'ISI'  ...
    'BurstDur'    'OnsetDelay'     'RiseDur'           'FallDur'           'ITD'   ...
    'ITDtype'    'TimeWarpFactor'  'Fsam'              'Fcar'              'Fmod' ...
    'LowCutoff'    'HighCutoff'    'FreqWarpFactor'     'CarStartPhase'    'ModStartPhase'    ...
    'ModDepth'     'SPL'           'SPLtype'            'DAC'              'Ntone' ...
    'Nrep'         'ModTheta'};


Missing = setdiff(RFNS, fieldnames(p));
FixStr = [char(10) '========> stimparams' p.StimType ' needs to be fixed (see GenericStimparams).'];
if ~isempty(Missing),
    Mess = ['Missing stim parameter(s): ' cell2words(Missing) '.' FixStr];
    return;
end
% test for sizes and types
i_scalar = [2 3 4 13 25 26];
i_char = [1 11 23 24];
i_3D = [14 19];
i_dual = [6 8 9 15 16 17 20 21 22 27];
i_single = setdiff(1:numel(RFNS), [i_scalar i_char i_3D i_dual]);
Ncond = p.Ncond;
Ntone = p.Ntone;
for ii=i_scalar,
    fn = RFNS{ii};
    V = p.(fn);
    if ~isnumeric(V) || ~isequal(1,numel(V)),
        Mess = ['Generic stimulus parameter ''' fn ''' must be single number.' FixStr];
        return;
    end
end
for ii=i_char,
    fn = RFNS{ii};
    V = p.(fn);
    if ~ischar(V) || ~isequal(1,size(V,1)),
        Mess = ['Generic stimulus parameter ''' fn ''' must be horizontal char string.' FixStr];
        return;
    end
end
for ii=i_dual,
    fn = RFNS{ii};
    V = p.(fn);
    if ~isnumeric(V) || ~isequal([Ncond 2],size(V)) 
        Mess = ['Generic stimulus parameter ''' fn ''' must be (Ncond x 2) sized numeric array.' FixStr];
        return;
    end
end
for ii=i_single,
    fn = RFNS{ii};
    V = p.(fn);
    if ~isnumeric(V) || ~isequal([Ncond 1],size(V))
        Mess = ['Generic stimulus parameter ''' fn ''' must be (Ncond x 1) sized numeric array.' FixStr];
        return;
    end
end
for ii=i_3D,
    fn = RFNS{ii};
    V = p.(fn);
    Sz = size(V); 
    if isequal(1,Ntone), Sz = [Sz 1]; end
    if ~isnumeric(V),
        Mess = ['Generic stimulus parameter ''' fn ''' must be numeric array.' FixStr];
        return;
    end
    if ~isequal([Ncond 2 Ntone], Sz),
        Mess = ['Generic stimulus parameter ''' fn ''' must have size [Ncond x 2] or [Ncond x 2 x Ntone].' FixStr];
        return;
    end
end
% special case; OnsetDelay. Must correspond to integer # samples
dt = 1e3/p.Fsam;
N = p.OnsetDelay/dt; 
N = N(N~=0); % restrict test to nonzero delays
RD = abs(reladev(N,round(N)));
if any(RD>1e-6),
    Mess = ['Generic stimulus parameter ''OnsetDelay'' must correspond to an integer number of samples.' FixStr];
    return;
end
% all is well if we get here
P = p;



