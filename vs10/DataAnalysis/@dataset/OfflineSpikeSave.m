function Mess = OfflineSpikeSave(DS, SPT, DetectCall, AddInfo);
% Dataset/OfflineSpikeSave - save spike arrival times from offline analysis
%    Mess = OfflineSpikeSave(D, SPT, DetectCall, S) saves spiketimes SPT 
%    extracted from dataset D. In case of a problem, Mess contains a
%    message identifying the problem ('' when okay).
%    INPUTS
%         D: dataset from which the spiketimes were extracted.
%       SPT: Ncond x Nrep cell array. SPT{icond, irep} is a row array 
%            holding the spike arrival times of condition icond, repetition 
%            irep, in ms after stimulus onset.
%   DetectCall: cell array allowing the recalculation of the detection.
%            This cell array has the form {@Foo Dum Arg1 Arg2 ...}, where
%               @Foo is a function handle to the spike-detecting function
%                Dum is an arbitrary dummy value (a stand-in for D)
%                Arg1 ... are additional input arguments to Foo.
%             Calling @Foo(D, Arg1, ...) should result in the recalculation
%             of spike detection, and Foo is supposed to save the spikes
%             by calling OfflineSpikeSave.
%         S: a struct containing any information you want to save along with
%            the actual spike times. It can be retrieved by
%            Dataset/OfflineSpikeInfo.
%    The spike times are saved in the 'Early_Spiketimes' subdir of 
%    processed_dir(). This folder must exist prior to the call. Within it,
%    each experiment has its own subdir which will be created
%    automatically. Within the subdirs, each dataset has its own file. 
%    In an obvious notation, the full filename is 
%
%      <processed_datadir>\Early_spiketimes\<Expname>\DSxxxx.offspt
%    
%    An error is thrown when saving of the spike times fails.
%
%    Example code
%      function MySpikeDetector(D, X, Y, Z) % D is dataset; X,Y,Z detection parameters
%      ..... % here the spikes are detected and placed in cell array SPT
%      S.Note = 'Reliable spike detection.';
%      OfflineSpikeSave(D, SPT, {fhandle(mfilename), [], X, Y, Z}, S.Note)
%
%    See also Dataset/OfflineSpikeRemove, Dataset/HasOfflineSpikes,
%    Dataset/OfflineSpikeRecalc, Dataset/OfflineSpikeList,
%    Dataset/OfflineSpikeInfo.

if ~isstruct(AddInfo),
    error('Fourth argument of OfflineSpikeSave should be struct.');
end

% Check type & size and content of SPT
szSPT = size(SPT);
NcNr = [Ncond(DS), repcount(DS)];
isrow = @(x)(isvector(x)&&size(x,1)<2)||isempty(x);
anisnan = @(x)any(isnan(x));
if ~iscell(SPT),
    error('SPT input argument must be a Ncond x RepCount cell array.');
elseif ~isequal(NcNr, szSPT),
    error(['Incorrect size of spiketime cell array; should be Ncond x RepCount = ' sizeString(NcNr) '.']);
elseif ~all(all(cellfun(@isnumeric, SPT))),
    error('All cells of SPT must be numeric.');
elseif ~all(all(cellfun(isrow, SPT))),
    error('Each cell of SPT must be either a row vector or empty.');
elseif any(any(cellfun(anisnan, SPT))),
    error('Spiketimes may not contain NaNs.');
end
% check type and content of P
if ~iscell(DetectCall),
    error('Input arg DetectCall must be cell array.');
elseif ~isfhandle(DetectCall{1}),
    error('First element of cell array DetectCall must be a function handle.');
end

% ID info
Expname = expname(DS);
iRec = irec(DS);

% save
ID = CollectInStruct(Expname, iRec);
S = structJoin(ID, CollectInStruct('-Recall', DetectCall), CollectInStruct('-SpikeTimes', SPT), ...
    CollectInStruct('-Additional_Info', AddInfo));
FN = OfflineSpikeFilename(DS,true); % 2nd arg: create exp subdir if needed
Mess = '';
try,
    structSave(FN, S, '-V6'); % version 6: disable time consuming data compression - not worth it.
catch,
    Mess = ['Problem while saving spike times.' char(10) lasterr];
end


















