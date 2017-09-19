function D = ZWOAEimport(ExpID, idata, flag)
% ZWOAEimport - retrieve recorded ZWOAE data from file
%
%   D = ZWOAEimport is a shell around the function getZWOAEdata. It uses
%   getZWOAEdata to actually read the data, but differs from getZWOAEdata
%   in the following ways:
%       - the time signal is replaced by a spectrum
%       - caching is used
%       - fieldnames are changed for better overview of the data.
%
%   When the recording was 'CM' (see field D.RecType) the amplitude spectrum
%   is in dB re 1V, rather than in dB SPL (see also correctCMamplitude)
%
%     D = ZWOAEimport(47, 154) import data from gerbil #47, recording 154 (see getZWOAEdata)
%     D = ZWOAEimport(47, [1 2 4]) returns datasets #1, 2 and 4 in a struct array.
%
%     D = ZWOAEimport(.. , '-nosig') removes the spectrum from the datastruct
%       and just returns the measurement parameters; this is useful for
%       bookkeeping purposes.
%
%     D = ZWOAEimport(.. , '-keeptimesig') keeps the (loop-averaged) time
%     signal in addition to its spectrum.
%
%  See also getZWOAEdata, ZWOAEspec, ZWOAEfilename.

if nargin<3, flag=''; 
elseif ~isempty(flag),
    [flag, Mess] = keywordMatch(flag, {'-nosig', '-keeptimesig'}, 'flag argument');
    error(Mess);
end

if ~isnumeric(ExpID),
    error('First input arg must be gerbil number, not exp name.');
end

if numel(idata)>1, % multiple datasets: recursion
    for ii=1:numel(idata);
        D(ii) = ZWOAEimport(ExpID, idata(ii), flag);
    end
    return;
end

%=========single file from here=============

%-------imaginary dataset?
if ~isreal(idata),
    D = imagZWOAErecording(ExpID, idata);
    if isequal('-nosig', flag),
        D = rmfield(D, {'sp__________________', 'df', 'MG', 'PH', 'Nsam', 'Pnf'});
    elseif isequal('-keeptimesig', flag),
        error('Invalid use of -keeptimesig flag with pooled data.');
    end
    return;
end

% caching (only on single datasets!)
CFN = [mfilename '_' int2str(ExpID) flag];
if isequal('-nosig', flag), Ncache=1e4;
else, Ncache = 500;
end
CacheParam = {ExpID, idata};
D = FromCacheFile(CFN,CacheParam);
if ~isempty(D), return; end % cached, we're done

% ----read data and shape them----------
DD = getZWOAEdata(ExpID, idata);

% Two options, either ExpID is an 'old' recording (pre gerbil 74)
% or a new one. The new ones are recognized by their fieldnames
% Old and new recordings have to be shaped differently to get the same
% output format. The new ones are handled by a local function, the old ones
% are "as is"

if isfield(DD, 'Stimparam'),
    %it is a new one, delegate to local function.
    D = local_new2D(DD, flag);
else,  %old one (igerbil<74), proceed as always
    %=====GENERIC RECORDING PARAMS=======
    fs = DD.fs;
    periodicity = DD.periodicity;
    micgain = DD.micgain;
    dBV = DD.micsensitivity.dBV;
    dBSPL = DD.micsensitivity.dBSPL;
    %=====SPECIFIC RECORDING PARAMS=======
    recID = DD.index; 
    RecType = getFieldOrDefault(DD.setupinfo,'rectype','unknown');
    companionID = getFieldOrDefault(DD.setupinfo,'companionID',[]);
%    companionID = DD.setupinfo.companionID;
    % =========STIMULUS PARAMS=========
    %==Levels
    L1 = DD.stimpars.L1;
    L2 = DD.stimpars.L2;
    Lsup = getFieldOrDefault(DD.stimpars,'Lsup',-999); 
    Lmim = -999; 
    %==Frequencies, still in Hz
    F1 = DD.stimpars.F1; % Hz
    F2 = DD.stimpars.F2; % Hz
    Fsup = getFieldOrDefault(DD.stimpars,'Fsup',999); % Hz
    Fmim = 0; 
    %==Phases; still in rad
    PH1 = getFieldOrDefault(DD.stimpars,'ph1',999); % rad
    PH2 = getFieldOrDefault(DD.stimpars,'ph2',999); % rad
    PHsup = getFieldOrDefault(DD.stimpars,'phsup',-999);  % rad
    PHmim = 999;
    %Frequency unit conversion: Hz --> kHz
    if ~isequal(F1, 999), F1 = F1/1e3; end
    if ~isequal(F2, 999), F2 = F2/1e3; end
    if ~isequal(Fsup, 999), Fsup = Fsup/1e3; end
    %Phase unit conversion: rad --> cycle
    if ~isequal(PH1, 999), PH1 = PH1/2/pi; end
    if ~isequal(PH2, 999), PH2 = PH2/2/pi; end
    if ~isequal(PHsup, 999), PHsup = PHsup/2/pi; end
    %give probename a value that tells software to ignore it
    probename = 'none';
    
    % =====FIGURE OUT WHETHER F1 or F2 is ZWUIS=====
    num1 = numel(F1); num2 = numel(F2);
    if num1 > 1, % F1 is ZWUIS
        if num2 > 2, % F2 is also ZWUIS --> error!
            error('F1 and F2 are both ZWUIS, not yet implemented');
        end
        % F1 = ZWUIS, F2 = single
        Fzwuis = F1; Lzwuis = L1; PHzwuis = PH1;
        Fsingle = F2; Lsingle = L2; PHsingle = PH2;
    else, % F2 is ZWUIS, F1 = single
                   % or neither is ZWUIS
        Fzwuis = F2; Lzwuis = L2; PHzwuis = PH2;
        Fsingle = F1; Lsingle = L1; PHsingle = PH1;
    end
    Nzwuis = numel(Fzwuis); % # of ZWUIS components
    Fz_mean = mean(Fzwuis); % mean ZWUIS frequency

    %---determine whether ZWUIS is 'above','below' or 'none'---
    ZwuisLoc = 'above';
    if Nzwuis == 1;  % 2-tone
        ZwuisLoc = 'none';
    elseif Fsingle > Fz_mean  % ZWUIS below
        ZwuisLoc = 'below';
    end
    %======CLEAN UP THE STIMULUS PARAMETERS=========
    %---discard "soft" components by assigning approp. values---
    %NOTE: numel(Lzwuis), numel(Lsingle), numel(Lsup) == 1 
    if Lzwuis < -100,
        Lzwuis = -999; % set level to -999
        Fzwuis = 0; % set frequency to 0
        PHzwuis = 999; % set phase to 999
    end
    if Lsingle < -100,
        Lsingle = -999;
        Fsingle = 0;
        PHsingle = 999;
    end
    if Lsup < -100,
        Lsup = -999;
        Fsup = 0;
        PHsup = 999;
    end

    signal = DD.signal; %get averaged signal
    
    %==================================================
    %==================================================
    %all params have been found, collect in structure D
    D = CollectInStruct(ExpID, recID, RecType, companionID, ZwuisLoc, Nzwuis, fs, periodicity, '-', Fzwuis,...
        Fz_mean, Fsingle, Fsup, Fmim, Lzwuis, Lsingle, Lsup, Lmim, PHzwuis,...
        PHsingle, PHsup, PHmim, '-',micgain, dBV, dBSPL, probename, signal);
    
end %this is the end of processing of either new (post gerbil 73) or old (pre gerbil 74) recording

D.Nzwuis = numel(D.Fzwuis);

% what to do with the signal depends on flag
if isequal('', flag) || isequal('-keeptimesig', flag), % compute spectrum
    D.sp__________________ = '________________';
    [D.df, D.MG, D.PH, D.Nsam, D.Pnf] = ZWOAEspec(D); %
    %correction of CM amplitudes
    try,
        D = correctCMamplitude(D);
    catch,
        disp('CM magnitude spectrum was not corrected upon import.');
    end
end
% unless told to keep time signal, toss it
if ~isequal('-keeptimesig',flag), D = rmfield(D,'signal'); end

%anecdotal corrections
D = local_correctJunk(D);

ToCacheFile(CFN, Ncache, CacheParam, D);

%====LOCAL===

function D = local_new2D(DD, flag)

%get the units specified in the saved data
unitF = DD.Userparam.unitF;
unitL = DD.Userparam.unitL;
unitPH = DD.Userparam.unitPH;
%use them to calculate conversion factors:
%Freqs must be in kHz
%Levels must be in dB SPL
%Phases must be in cycles
[Ffac,Mess] = UnitConvert(unitF, 'kHz'); error(Mess);
[Pfac, Mess] = UnitConvert(unitPH, 'cycle');  error(Mess);
if ~isequal('dB SPL', unitL),
    fn = ZWOAEfilename(DD.Expinfo.ExpID, DD.Other.recID);
    error([fn ' : recorded levels are not in dB SPL. This is required.']);
end

%---General info---
ExpID = DD.Expinfo.expID;
recID = DD.Other.recID;
RecType = DD.Other.recType;
companionID = DD.Other.companionID;
Nzwuis = DD.Stimparam.Nzwuis;
periodicity = DD.Stimparam.avsamp;
%---Frequencies (including sample frequency), must be in kHz---
fs = DD.Stimparam.Fsam*Ffac;
Fzwuis = DD.Stimparam.Fzwuis*Ffac;
Fz_mean = DD.Stimparam.Fz_mean*Ffac; 
Fsingle = DD.Stimparam.Fsingle*Ffac;
Fsup = DD.Stimparam.Fsup*Ffac;
Fmim = getFieldOrDefault(DD.Stimparam, 'Fmim', 0);
%---Levels, must be in dB SPL---
Lzwuis = DD.Stimparam.Lzwuis;
Lsingle = DD.Stimparam.Lsingle;
Lsup = DD.Stimparam.Lsup;
Lmim = getFieldOrDefault(DD.Stimparam, 'Lmim', -999);
%---PHase, must be in cycle---
PHzwuis = DD.Stimparam.PHzwuis*Pfac;
PHsingle = DD.Stimparam.PHsingle*Pfac;
PHsup = DD.Stimparam.PHsup*Pfac;
PHmim = getFieldOrDefault(DD.Stimparam, 'PHmim', 999);
%---Microphone info---
micgain = DD.Setupinfo.micgain;
dBV = DD.Setupinfo.dBV;
dBSPL = DD.Setupinfo.dBSPL;

%---probe name for inclusion of its transfer function
probename = getFieldOrDefault(DD.Setupinfo, 'probename', 'none');

%Find out whether ZWUIS is below, above or none (e.g. 2-tone stimulus)
ZwuisLoc = 'above';
if Nzwuis == 1;  %2-tone
    ZwuisLoc = 'none';
elseif Fsingle > Fz_mean  %ZWUIS below
    ZwuisLoc = 'below';
end

%get the averaged signal, needed for calculation of spectrum
signal = DD.Data.processed;

%collect everything in struct D
D = CollectInStruct(ExpID, recID, RecType, companionID, ZwuisLoc, Nzwuis, fs, periodicity, '-', Fzwuis,...
    Fz_mean, Fsingle, Fsup, Fmim, Lzwuis, Lsingle, Lsup, Lmim, PHzwuis,...
    PHsingle, PHsup, PHmim, '-',micgain, dBV, dBSPL, signal,probename);

function D = local_correctJunk(D);
% correct incidental bookkeeping errors
if isequal(80,D.ExpID) & ismember(D.recID,2230:2242), % these are Acoustic recordings, not CM 
    D.RecType = 'ACOUSTIC';
end



