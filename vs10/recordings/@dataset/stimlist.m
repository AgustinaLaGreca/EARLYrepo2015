function [L,R, T] = stimlist(D, flag)
% Dataset/IDstring - stimulus info in a struct.
%   stimlist(DS) returns a struct with char-valued fields
%         iRec: # recording in experiment
%     IDstring: identifier returned by IDstring(DS)
%       Xrange: range of X variable
%          SPL: intensities
%        Spec1: special parameter #1
%        Spec2: special parameter #2
%          Dur: string like '12x10x70/100 ms'
%         Chan: active DA channels
%
%   [L, R, T] = stimlist(DS) also returns recording info in struct array R 
%   and stimulus type in cell string T. Dataset arrays produce struct array
%   output.
%
%   stimlist(DS,'clear') clears the cache for the experiment to which DS
%   belongs.
%      
%   Stimlist is useful for listing data of an experiment. It is used by
%   databrowse.
%
%   See also Dataset/paramview, Dataset/IDstring, Experiment/stimparam.

flag = arginDefaults('flag', '');

CFN = [mfilename '_' expname(D(1))]; % cache file name
if isequal('clear', flag),
    rmcache(CFN);
    return;
end

if numel(D)>1,
    for ii=1:numel(D),
        [L(ii), R(ii) T(ii)] = stimlist(D(ii));
    end
    L = reshape(L, size(D));
    R = reshape(R, size(D));
    T = reshape(T, size(D));
    return;
end

% single DS from here
STR = stimlist_strfun(D); % helper functions for string conversion
[LRT, CFN, CP] = getcache(CFN, irec(D));
if ~isempty(LRT), 
    [L,R, T] = deal(LRT{:});
    return; 
end

if isvoid(D),
    error('Void dataset has no stimulus.');
elseif isdummy(D,'nostim'), % retrieve complete version
    D = read(dataset(), expname(D), irec(D));
end

L.iRec = num2str(irec(D));
L.IDstring = IDstring(D);
L.Xrange = STR.xrange(D.Stim.Presentation.X);
L2 = struct('SPL', '', 'Spec1', '', 'Spec2', '', 'Dur', '', 'Chan', '');
switch lower(stimtype(D)),
    case 'binzw',
        L2.SPL = STR.shstring(D.Stim.SPL); % level per component
        L2.SPL = [L2.SPL ' dB/cmp'];
        L2.Spec1 = [STR.shstring(D.Stim.LowFreq) '-' STR.shstring(D.Stim.HighFreq) ' Hz'];  % zwuis freq
        L2.Spec2 = [STR.shstring(D.Stim.Ncomp) ' cmp' ];  % zwuis freq
    case 'fs',
        L2.SPL = [STR.shstring(D.Stim.SPL) ' dB']; % tone level
        L2.Spec1 = STR.modstr(D.Stim); % mod
    case 'mask',
        L2.SPL = [STR.shstring(D.Stim.MnoiseSPL) ' dB']; % noise level
        L2.Spec1 = [STR.shstring(D.Stim.ToneFreq) '-Hz tone']; % tone freq
        L2.Spec2 = ['S/N=' STR.shstring(D.Stim.ToneSPL) ' dB']; % S/N
    case {'nphi' 'nphi2'},
        L2.SPL = [STR.shstring(D.Stim.SPL) ' ' D.Stim.SPLUnit]; % noise level
        L2.Spec1 = STR.modstr(D.Stim); % mod
        L2.Spec2 = [STR.shstring([D.Stim.LowFreq D.Stim.HighFreq], '..') ' Hz']; % noise cutoffs
    case 'rf',
        L2.SPL = [STR.shstring(D.Stim.StartSPL) ':' STR.shstring(D.Stim.StepSPL) ':' STR.shstring(D.Stim.EndSPL) 'dB SPL']; % probe level
        L2.Spec1 = STR.modstr(D.Stim); % mod
    case 'sup',
        L2.SPL = [STR.shstring(D.Stim.SPL) '-dB probe']; % probe level
        L2.Spec1 = [STR.shstring(D.Stim.SupCenterFreq) '-Hz, ' num2str(D.Stim.SupNcomp) '-comp'];  % suppr freq
        L2.Spec2 = [STR.shstring(D.Stim.SupSPL) '-dB sup']; % suppr level
    case 'irn'
        L2.SPL = [STR.shstring(D.Stim.SPL) 'dB']; % probe level
        L2.Spec1 = [STR.shstring([D.Stim.Startdelta_T D.Stim.Enddelta_T],'..') 'ms'];  % suppr freq
        L2.Spec2 = [STR.shstring([D.Stim.LowFreq D.Stim.HighFreq], '..') ' Hz']; % suppr level 
    case 'hp'
        L2.SPL = [STR.shstring(D.Stim.SPL) 'dB']; % probe level
        L2.Spec1 = [STR.shstring([D.Stim.StartFreq D.Stim.EndFreq],'..') 'ms'];  % suppr freq
        L2.Spec2 = [STR.shstring([D.Stim.LowFreq D.Stim.HighFreq], '..') ' Hz']; % suppr level 
    case 'spont'
        L2.SPL = [STR.shstring(-Inf) 'dB']; % probe level
        L2.Spec1 = [STR.shstring([D.Stim.ISI]) 'ms'];  % suppr freq
        L2.Spec2 = ''; % suppr level 
    case 'zw',
        L2.SPL = STR.shstring(D.Stim.SPL); % level per component
        if ~isequal('-', D.Stim.ProfileType), % % also display SPL jump
            jump = D.Stim.ProfileSPLjump;
            if jump<0, L2.SPL = [L2.SPL  '-'  STR.shstring(abs(jump))];
            else, L2.SPL = [L2.SPL '+'  STR.shstring(abs(jump))];
            end
        end
        switch D.Stim.ProfileType,
            case 'single cmp', apx = ' 1';
            case 'low end', apx = ' >';
            case 'high end', apx = ' <';
            otherwise, apx = '';
        end
        L2.SPL = [L2.SPL ' dB/cmp ' apx];
        L2.Spec1 = [STR.shstring(D.Stim.LowFreq) '-' STR.shstring(D.Stim.HighFreq) ' Hz'];  % zwuis freq
        L2.Spec2 = [STR.shstring(D.Stim.Ncomp) ' cmp' ];  % zwuis freq
    otherwise, % more "object-oriented" approach: delegate to stimulus-specific function
        stimdisplayer = ['stimdisp' upper(stimtype(D))];
        if exist(stimdisplayer, 'file'),
            [L2.SPL L2.Spec1, L2.Spec2] = feval(stimdisplayer, D.Stim);
        else,
            warning(['No stimulus display function ''' stimdisplayer ''' found.'])
        end
end
L2.Dur = durationstring(D, 'compact');
L2.Chan = upper(D.Stim.DAC(1));
L = structJoin(L, L2);
R = local_R(D);
T = {stimtype(D)};
putcache(CFN, 1e3, CP, {L R T});

%=====================================
%=====================================
%=====================================
function Str = local_xr(D);
X = D.Stim.Presentation.X;
[Xmin, Xmax] = minmax(X.PlotVal);
U = X.ParUnit;
if isequal('octave', lower(U)),
    U = 'Oct';
end
if isequal('Hz', U),
    Str = [STR.num2kstr(Xmin) ':' STR.num2kstr(Xmax) ' ' U];
else,
    Str = [STR.shstring(Xmin) ':' STR.shstring(Xmax) ' ' U];
end

function StimStr = local_dur(D);
Dur = burstdur(D,1);
Dur = unique(0.1*round(10*mean(Dur,1)));
DurStr = [STR.shstring(Dur) '/' num2str(unique(0.1*round(10*mean(repdur(D),1)))) ' ms'];
Pres = D.Stim.Presentation;
if has2varparams(D),
    StimStr = [num2str(D.Stim.Ncond_XY(1)) 'x' num2str(D.Stim.Ncond_XY(2)) 'x' num2str(Pres.Nrep) 'x' DurStr];
else,
    StimStr = [num2str(Pres.Ncond) 'x' num2str(Pres.Nrep) 'x' DurStr];
end
if ~isempty(D.Stopped),
    StimStr = [StimStr ' ***'];
end


function R = local_R(D);
% assuming RX6 standard layout
if isfield(D.Data, 'RX6_digital_1'),
    R.DigRec = 'Dig-1';
else,
    R.DigRec = ' --  ';
end
A1 = isfield(D.Data, 'RX6_analog_1');
A2 = isfield(D.Data, 'RX6_analog_2');
if A1 && ~A2,
    Astr = 'AD-1  ';
elseif ~A1 && A2,
    Astr = 'AD-2  ';
elseif A1 && A2,
    Astr = 'AD-1&2';
else,
    Astr = '----  ';
end
R.AnaRec = Astr;






