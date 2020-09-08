function ArgOut = EvalSCCXCC(varargin)
%EVALSCCXCC  calculate the SCC, XCC and DIFCOR.
%   T = EVALSCCXCC(ds1P, SubSeq1P, ds1N, SubSeq1N, ds2P, SubSeq2P, ds2N, SubSeq2N)
%   calculates the SCC, XCC and DIFCOR from the specified responses. The noise token
%   A+ response of the first cell is given by subsequence SubSeq1P from dataset ds1P.
%   Response to the A- token is given by subsequence SubSeq1N of dataset ds1N. Mutatis
%   mutandis for the second cell.
%   E.g.
%           ds1p = dataset('A0241', -128); %75-3-NSPL-A+
%           ds1n = dataset('A0241', -129); %75-4-NSPL-A-
%           ds2p = dataset('A0241', -133); %76-3-NSPL-A+
%           ds2n = dataset('A0241', -134); %76-4-NSPL-A-
%           EvalSCCXCC(ds1p, 50, ds1n, 50, ds2p, 50, ds2n, 50);
%
%   T = EVALSCCXCC(ds1, SubSeqs1, ds2, SubSeqs2) calculates the SCC, XCC and DIFCOR 
%   where the responses of a cell to the two noise tokens is located in the same
%   dataset at different subsequence numbers.
%   E.g.
%           ds1 = dataset('A0242', 60); %15-2-NRHO-A+-
%           ds2 = dataset('A0242', 50); %13-2-NRHO-A+-
%           EvalSCCXCC(ds1, [+1, -1], ds2, [+1, -1]);
% 
%   Optional properties and their values can be given as a comma-separated
%   list. To view list of all possible properties and their default value, 
%   use 'factory' as only input argument.

%B. Van de Sande 18-04-2005

%Attention! Caching system is not supported because the MEX implementation of
%SPTCORR makes the overhead of caching unrewarding ...


%% ---------------- CHANGELOG -----------------------
%  Wed Apr 25 2012  Abel   
%   - Added support for gabor fitfnc
%  Fri May 29 2020 Darina
%   - Corrected the parsing of spike times cells

%-----------------------------------template---------------------------------
%Attention! The cochlear distance, which can be derived from the tuning frequency
%of a cell (either the characteristic frequency CF or the dominant frequency DF),
%is not included in the table. Using the function GREENWOOD.M, delivered with the
%SGSR distribution, this is easily calculated afterwards when using STRUCTTOOLS.
%For the same reason the difference in tuning frequency between the two crosscorrelated
%cells is not included ...
%... identification ...
Template.ds1p.filename       = '';        %Datafile name for datasets of first input
Template.ds1p.icell          = NaN;       %Cell number of first input
Template.ds1p.iseq           = NaN;       %Sequence number of first input's positive response 
Template.ds1p.seqid          = '';        %Identifier of first dataset
Template.ds1p.isubseq        = NaN;       %Subsequence number of spiketrain used for first dataset
Template.ds1n.iseq           = NaN;       %Sequence number of second dataset
Template.ds1n.seqid          = '';        %Identifier of second dataset
Template.ds1n.isubseq        = NaN;       %Subsequence number of spiketrain used for second dataset
Template.ds2p.filename       = '';        %Datafile name for second dataset
Template.ds2p.icell          = NaN;       %Cell number of second dataset
Template.ds2p.iseq           = NaN;       %Sequence number of second dataset
Template.ds2p.seqid          = '';        %Identifier of second dataset
Template.ds2p.isubseq        = NaN;       %Subsequence number of spiketrain used for second dataset
Template.ds2n.iseq           = NaN;       %Sequence number of second dataset
Template.ds2n.seqid          = '';        %Identifier of second dataset
Template.ds2n.isubseq        = NaN;       %Subsequence number of spiketrain used for second dataset
%... miscellaneous ...
Template.tag                 = 0;         %General purpose tag field
Template.createdby           = mfilename; %Name of MATLAB function that generated the data
%... stimulus parameters ...
%Stimulus parameters are saved as numerical matrices where different columns correspond to
%different datasets and different rows designate different channels ... 
Template.stim.burstdur       = NaN(1, 4); %Stimulus duration in ms
Template.stim.repdur         = NaN(1, 4); %Repetition duration in ms
Template.stim.nrep           = NaN(1, 4); %Number of repetitions
Template.stim.spl            = NaN(2, 4); %Sound pressure level in dB
Template.stim.avgspl1        = NaN;               %Averaging of SPL in the power spectrum for one
Template.stim.avgspl2        = NaN;               %input ...
Template.stim.avgspl         = NaN;               %Averaging of SPL in the power spectrum for both
                                                  %inputs ...
%... threshold curve ...
Template.thr1.cf             = NaN;       %Characteristic frequency retrieved from threshold curve
Template.thr1.sr             = NaN;       %Spontaneous rate retrieved from threshold curve
Template.thr1.thr            = NaN;       %Threshold at characteristic frequency
Template.thr1.q10            = NaN;       %Q10 retrieved from threshold curve
Template.thr1.bw             = NaN;       %Width 10dB above threshold (Hz)
Template.thr2.cf             = NaN;
Template.thr2.sr             = NaN;
Template.thr2.thr            = NaN;
Template.thr2.q10            = NaN;
Template.thr2.bw             = NaN;
%... rate curve ...
%Attention! Driven rate can be calculated afterwards by subtracting spontaneous rate derived
%from threshold curves from the mean rate ...
Template.rc1.mean            = NaN;
Template.rc2.mean            = NaN;
%... autocorrelogram ...
Template.sac1.max            = NaN;       %Maximum of shuffled autocorrelogram (DriesNorm)
Template.sac1.hhw            = NaN;       %Half height width on autocorrelogram (ms)
Template.sac1.fft.df         = NaN;       %Dominant frequency in autocorrelogram (Hz)
Template.sac1.fft.bw         = NaN;       %Bandwidth (Hz)
Template.dac1.max            = NaN;       %Maximum of diffautocorrelogram (DriesNorm)
Template.dac1.maxsecpks      = [NaN, NaN];%Height of secondary peaks (DriesNorm)
Template.dac1.lagsecpks      = [NaN, NaN];%Lag at secondary peaks (ms)
Template.dac1.fft.df         = NaN;       %Dominant frequency in diffcorrelogram (Hz)
Template.dac1.fft.bw         = NaN;       %Bandwidth (Hz)
Template.dac1.env.hhw        = NaN;       %Half height with on envelope of sumcorrelogram (ms)
Template.mac1.max            = NaN;       %Maximum of sumautocorrelogram (DriesNorm)
Template.mac1.maxsecpks      = [NaN, NaN];%Height of secondary peaks (DriesNorm)
Template.mac1.lagsecpks      = [NaN, NaN];%Lag at secondary peaks (ms)
Template.mac1.fft.df         = NaN;       %Dominant frequency in sumcorrelogram (Hz)
Template.mac1.fft.bw         = NaN;       %Bandwidth (Hz)
Template.mac1.env.hhw        = NaN;       %Half height with on envelope of sumcorrelogram (ms)
Template.sac2.max            = NaN;       
Template.sac2.hhw            = NaN;       
Template.sac2.fft.df         = NaN;      
Template.sac2.fft.bw         = NaN;       
Template.dac2.max            = NaN;       
Template.dac2.maxsecpks      = [NaN, NaN];
Template.dac2.lagsecpks      = [NaN, NaN];
Template.dac2.fft.df         = NaN;       
Template.dac2.fft.bw         = NaN;       
Template.dac2.env.hhw        = NaN;  
Template.mac2.max            = NaN;       
Template.mac2.maxsecpks      = [NaN, NaN];
Template.mac2.lagsecpks      = [NaN, NaN];
Template.mac2.fft.df         = NaN;       
Template.mac2.fft.bw         = NaN;       
Template.mac2.env.hhw        = NaN; 
%... crosscorrelogram ...
Template.scc.max             = NaN;       %Maximum of shuffled crosscorrelogram (DriesNorm)
Template.scc.rate            = NaN;       %Maximum of shuffled crosscorrelogram (Rate)
Template.scc.lagatmax        = NaN;       %Lag at maximum of SCC (ms)
Template.scc.maxsecpks       = [NaN, NaN];%Height of secondary peaks (DriesNorm)
Template.scc.lagsecpks       = [NaN, NaN];%Lag at secondary peaks (ms)
Template.scc.hhw             = NaN;       %Half height width on SCC (ms)
Template.scc.fft.df          = NaN;       %Dominant frequency in SCC (Hz)
Template.scc.fft.bw          = NaN;       %Bandwidth (Hz)
Template.dcc.max             = NaN;       %Maximum of diffcorrelogram (DriesNorm)
Template.dcc.rate            = NaN;       %Maximum of diffcorrelogram (Rate)
Template.dcc.lagatmax        = NaN;       %Lag at maximum (in ms)
Template.dcc.maxsecpks       = [NaN, NaN];%Height of secondary peaks (DriesNorm)
Template.dcc.lagsecpks       = [NaN, NaN];%Lag at secondary peaks (ms)
Template.dcc.fft.df          = NaN;       %Dominant frequecny in diffcorrelogram (Hz)
Template.dcc.fft.bw          = NaN;       %Bandwidth (Hz)
Template.dcc.env.max         = NaN;       %Maximum of enveloppe of diffcorrelogram (DriesNorm)
Template.dcc.env.lagatmax    = NaN;       %Lag at maximum of enveloppe (ms)
Template.dcc.env.hhw         = NaN;       %Half height with on envelope of diffcorrelogram (ms)
Template.mcc.rate            = NaN;       %Maximum of sumcorrelogram (Rate)
Template.mcc.lagatmax        = NaN;       %Lag at maximum (in ms)
Template.mcc.maxsecpks       = [NaN, NaN];%Height of secondary peaks (DriesNorm)
Template.mcc.lagsecpks       = [NaN, NaN];%Lag at secondary peaks (ms)
Template.mcc.fft.df          = NaN;       %Dominant frequecny in sumcorrelogram (Hz)
Template.mcc.fft.bw          = NaN;       %Bandwidth (Hz)
Template.mcc.env.max         = NaN;       %Maximum of enveloppe of sumcorrelogram (DriesNorm)
Template.mcc.env.lagatmax    = NaN;       %Lag at maximum of enveloppe (ms)
Template.mcc.env.hhw         = NaN;       %Half height with on envelope of sumcorrelogram (ms)
%by: Abel, gabor fit (see EvalNITD)
Template.gabor.bestitd        = NaN;
Template.gabor.bestitdc       = NaN;
Template.gabor.secpeaks       = [NaN, NaN];
Template.gabor.hhw            = NaN;
Template.gabor.hhwc           = NaN;
Template.gabor.envmax         = NaN;
Template.gabor.envpeak        = NaN;
Template.gabor.envpeakc       = NaN;
Template.gabor.carrfreq       = NaN;
Template.gabor.accfrac        = NaN;

%-------------------------------default parameters---------------------------
%Syntax parameters ...
DefParam.subseqinput   = 'indepval'; %'indepval' or 'subseq' ...
%Calculation parameters ...
DefParam.anwin         = [0 +Inf];   %in ms (Infinite designates stimulus duration) ...
DefParam.corbinwidth   = 0.05;       %in ms ...
DefParam.cormaxlag     = 15;         %in ms ...   
DefParam.envrunavunit  = '#';        %'#' or 'ms' ...
DefParam.envrunav      = 1;          %in ms or number of periods ...
DefParam.corfftrunav   = 100;        %in Hz ...
DefParam.diffftrunav   = 100;        %in Hz ...
DefParam.sumfftrunav   = 100;        %in Hz ...
DefParam.calcdf        = NaN;        %in Hz, NaN (automatic), 'cf' or 'df' ...
%Calculation of average SACs is only useful when supplying the responses from the same
%fiber or cell to two different noise tokens (e.g. A+, A-, B+, B-). The average SAC and
%XAC is calculated, but the results are stored in the returned structure-array using
%the fieldnames appropriate for SCCs and XCCs analysis. The generated plot also has this
%naming abnormality ...
% E.g.:
%       ds = dataset('R99040', '6-2-ab');
%       EvalSCCXCC(ds, 1, ds, 3, ds, 2, ds, 4, 'calctype', 'avgsac');
DefParam.calctype      = 'scc';      %'scc' or 'avgsac' ...
%Plot parameters ...
DefParam.plot          = 'yes';      %'yes' or 'no' ...
DefParam.corxrange     = [-5 +5];    %in ms ...
DefParam.corxstep      = 1;          %in ms ...
DefParam.fftxrange     = [0 500];    %in Hz ...
DefParam.fftxstep      = 50;         %in Hz ...
DefParam.fftyunit      = 'dB';       %'dB' or 'P' ...
DefParam.fftyrange     = [-50 10]; 
%by Abel: gabor fit 
DefParam.gaborfit    = 'no';
DefParam.samplerate  = 0.5;           %in number of elements per microsecond ...
DefParam.fitrange    = [];            %[-5000 +5000]; %in microseconds ... By default we use corxrange
DefParam.calcdf      = NaN;           %in Hz, NaN (automatic), 'cf' or 'df' ...


%----------------------------------main program------------------------------
%Evaluate input arguments ...
if (nargin == 1) && ischar(varargin{1}) && strcmpi(varargin{1}, 'factory')
    if (nargout == 0)
        disp('Properties and their factory defaults:');
        disp(DefParam);
    else
        ArgOut = DefParam;
    end
    return;
elseif (nargin == 2) && ischar(varargin{1}) && strcmpi(varargin{1}, 'checkprops') && isstruct(varargin{2})
    CheckParam(varargin{2});
    return;
else
    [Spt, Info, StimParam, Param] = ParseArgs(DefParam, varargin{:});
end

%Retrieve and calculate threshold curve information ...
ds1p = read(dataset, Info.ds1p.filename, Info.ds1p.iseq);
Thr(1) = getThr4Cell(ds1p.Experiment, Info.ds1p.icell);
ds2p = read(dataset, Info.ds2p.filename,Info.ds2p.iseq);
Thr(2) = getThr4Cell(ds2p.Experiment, Info.ds2p.icell); 
%fast and ugly fix for compatibility
for n = 1:length(Thr)
    if iscell(Thr(n).str)
        Thr(n).str = Thr(n).str{:};
    end
end

%Extract rate curve information ...
RC(1) = CalcRC(Spt{1:2}, Param);
RC(2) = CalcRC(Spt{3:4}, Param);

%Calculate correlograms ...
[SXAC(1), DAC(1), MAC(1)] = CalcAC(Spt{1:2}, Thr(1), Param);
[SXAC(2), DAC(2), MAC(2)] = CalcAC(Spt{3:4}, Thr(2), Param);
if strncmpi(Param.calctype, 's', 1)
    [SXCC, DCC, MCC] = CalcCC(Spt{:}, Thr, Param);
else
    [SXCC, DCC, MCC] = CalcAvgAC(Spt{:}, Thr, Param);
end

%Fit gabor function on CROSS correlation (see EvalNITD)
GBOR = [];
Param.gaborfit = strcmpi(Param.gaborfit, 'yes') || strcmpi(Param.gaborfit, 'y');
if Param.gaborfit
	% Determine frequency
	% - check input ([freq], 'cf' or 'df')
	if ~(isnumeric(Param.calcdf) && ((Param.calcdf > 0) || isnan(Param.calcdf))) && ...
			~(ischar(Param.calcdf) && any(strcmpi(Param.calcdf, {'cf', 'df'})))
		error('Property calcdf must be positive integer, NaN, ''cf'' or ''df''.');
	end
	% - set freq
	if isnumeric(Param.calcdf)
		%if calcdf was not NaN, just take the users input
		if ~isnan(Param.calcdf)
			Freq = Param.calcdf;
		%if is NaN and we have both THR's, take their average
		elseif ~isempty(Thr) && ~isnan(Thr(1).cf) && ~isnan(Thr(2).cf)
			Freq = (Thr(1).cf + Thr(2).cf)/2;
		%if is NaN and no THR available,
		elseif ~isempty(DCC.fft) && ~isnan( DCC.fft.df)
			Freq = DCC.fft.df;
		else
			Freq = SXCC.fft.df;
		end
	elseif strcmpi(Param.calcdf, 'cf')
		if ~isempty(Thr) && ~isnan(Thr(:).cf)
			Freq = (Thr(1).cf + Thr(2).cf)/2;
		else
			Freq = NaN;
		end
	elseif strcmpi(Param.calcdf, 'df')
		if ~isempty(DCC.diff) && ~isnan(DCC.fft.df)
			Freq = DCC.fft.df;
		else
			Freq = SXCC.fft.df;
		end
	else
		Freq = NaN;
	end
	
	% - Calc gabor
    gaborFailed = false;
	gaborErrStr = [];
    %bug in isnan metlab2009a
    if ~isempty(find(isnan(DCC.normco), 1)) || isempty(DCC.normco)
        gaborErrStr = 'diffcorr not available';
        gaborFailed = true;
    else
        try 
            GBOR = calcGABOR(DCC,Param,Freq);
        catch gaborErr
            gaborErrStr = gaborErr.message;
            gaborFailed = true;
        end
    end
	if gaborFailed
        warning('GABOR Fit failed => %s', gaborErrStr);
	end
	GBOR.err = gaborErrStr;
end

%debug
% sprintf(
% 'THR info 1:',
% Thr(1)
% Info.ds1p.filename
% Info.ds1p.icell
% Thr(2)
% Info.ds2p.filename
% Info.ds2p.icell  

%Display data ...
if strcmpi(Param.plot, 'yes')
    PlotData(SXAC, DAC, MAC, SXCC, DCC, MCC, GBOR, Thr, RC, Info, StimParam, Param);
end

%Return output if requested ...
if (nargout > 0)
    CalcData = Info; CalcData.stim = StimParam;
    [CalcData.rc1,  CalcData.rc2]  = deal(RC(1), RC(2));
    [CalcData.thr1, CalcData.thr2] = deal(Thr(1), Thr(2));
    [CalcData.sac1, CalcData.sac2] = deal(SXAC(1), SXAC(2));
    [CalcData.dac1, CalcData.dac2] = deal(DAC(1), DAC(2));
    [CalcData.mac1, CalcData.mac2] = deal(MAC(1), MAC(2));
    [CalcData.scc,  CalcData.dcc, CalcData.mcc]  = deal(SXCC, DCC, MCC);
	if Param.gaborfit && ~gaborFailed
		CalcData.gabor = GBOR;
	end
    ArgOut = structtemplate(CalcData, Template);
    % Add stuff that wouldn't work with the template ...
    ArgOut.dac1.x = DAC(1).lag;
    ArgOut.dac1.y = DAC(1).normco;
    ArgOut.dac2.x = DAC(2).lag;
    ArgOut.dac2.y = DAC(2).normco;
    ArgOut.dcc.x = DCC.lag;
    ArgOut.dcc.y = DCC.normco;
    ArgOut.mac1.x = MAC(1).lag;
    ArgOut.mac1.y = MAC(1).normco;
    ArgOut.mac2.x = MAC(2).lag;
    ArgOut.mac2.y = MAC(2).normco;
    ArgOut.mcc.x = MCC.lag;
    ArgOut.mcc.y = MCC.normco;
end

%----------------------------------------------------------------------------
function [Spt, Info, StimParam, Param] = ParseArgs(DefParam, varargin)

%Checking input arguments ...
Nds = length(find(cellfun('isclass', varargin, 'dataset')));
if (Nds == 2) %T = EVALSCCXCC(ds1, SubSeqs1, ds2, SubSeqs2)
    if (length(varargin) < 4)
        error('Wrong number of input arguments.');
    end
    if ~all(cellfun('isclass', varargin([1, 3]), 'dataset'))
        error('First and third argument should be datasets.'); 
    else
        [ds1p, ds1n, ds2p, ds2n] = deal(varargin{[1, 1, 3, 3]});
    end
    
    if ~isnumeric(varargin{2}) || ~any(length(varargin{2}) == [1, 2]) || ...
            ~isnumeric(varargin{4}) || ~any(length(varargin{4}) == [1, 2])
        error('Second and fourth argument should be scalars or two-element numeric vectors.')
    else
        InputVec = cat(2, varargin{2}([1, end]), varargin{4}([1, end]));
        InputVec(find(cellfun('length', varargin([2, 4])) == 1)*2) = NaN;
    end
    
    if isnan(InputVec(2)), ds1n = dataset; end
    if isnan(InputVec(4)), ds2n = dataset; end
    
    ParamIdx = 5;
elseif (Nds == 4) %T = EVALSCCXCC(ds1P, SubSeq1P, ds1N, SubSeq1N, ds2P, SubSeq2P, ds2N, SubSeq2N)
    if (length(varargin) < 8), error('Wrong number of input arguments.'); end
    
    if ~all(cellfun('isclass', varargin(1:2:7), 'dataset') | cellfun('isclass', varargin(1:2:7), 'edfdataset'))
        error('First, third, fifth and seventh argument should be datasets.');
    else
        [ds1p, ds1n, ds2p, ds2n] = deal(varargin{1:2:7});
    end

    if ~all(cellfun('isclass', varargin(2:2:8), 'double')) || ~all(cellfun('length', varargin(2:2:8)) == 1)
        error('Subsequences should be specified using a numerical scalar.');
    else
        InputVec = cat(2, varargin{2:2:8});
    end
    
    if isvoid(ds1n) && ~isnan(InputVec(2))
        error('Associated subsequence of void dataset should be NaN.');
    end
    if isvoid(ds2n) && ~isnan(InputVec(4))
        error('Associated subsequence of void dataset should be NaN.');
    end
        
    if (~isvoid(ds1n) && (~isequal(ds1p.filename, ds1n.filename) || ~isequal(ds1p.icell, ds1n.icell))) || ...
       (~isvoid(ds2n) && (~isequal(ds2p.filename, ds2n.filename) || ~isequal(ds2p.icell, ds2n.icell)))
        error('Responses to different noise tokens from same cell should also be recorded from same cell.'); 
    end
    
    ParamIdx = 9;
else
    error('Invalid input arguments.');
end

%Retrieving properties and checking their values ...
Param = checkproplist(DefParam, varargin{ParamIdx:end});
CheckParam(Param);

%Checking subsequences numbers and values of independent variable ...
dsNames = {'ds1p', 'ds1n', 'ds2p', 'ds2n'}; Nds = 4;
if strcmpi(Param.subseqinput, 'subseq')
    iSubSeqs = InputVec; IndepVals = NaN*zeros(1, Nds);
    for n = find(~isnan(iSubSeqs)) 
       try IndepVals(n) = eval(sprintf('%s.indepval(iSubSeqs(%d));', dsNames{n}, n)); 
       catch, error('One of the supplied subsequence numbers is invalid'); end
    end
else
    IndepVals = InputVec; iSubSeqs = NaN*zeros(1, Nds);
    for n = find(~isnan(IndepVals))
       idx = eval(sprintf('find(%s.Stim.Presentation.X.PlotVal == IndepVals(%d));', dsNames{n}, n));
       if ~isempty(idx) && (length(idx) == 1)
           iSubSeqs(n) = idx;
       else
           error('One of the supplied values of the independent variabale doesn''t exist.');
       end
    end
end

%Assembling spiketrains ...
Spt = cell(Nds, 1);
for n = find(~isnan(iSubSeqs))
    Spt{n} = eval(sprintf('%s.spiketimes(iSubSeqs(%d), :);', dsNames{n}, n));
end

%Assembling dataset information ...
Info.ds1p.filename  = lower(ds1p.ID.Experiment.ID.Name);
Info.ds1p.icell     = ds1p.ID.iCell;
Info.ds1p.iseq      = ds1p.ID.iDataset;
Info.ds1p.seqid     = lower([num2str(ds1p.ID.iCell) '-' num2str(ds1p.ID.iRecOfCell) '-' ds1p.StimType]);
Info.ds1p.isubseq   = iSubSeqs(1);
Info.ds1p.indepval  = IndepVals(1);
Info.ds1p.indepunit = ds1p.Stim.Presentation.X.ParUnit;

if ~isvoid(ds1n)
    Info.ds1n.iseq      = ds1n.ID.iDataset;
    Info.ds1n.seqid     = lower([num2str(ds1n.ID.iCell) '-' num2str(ds1n.ID.iRecOfCell) '-' ds1n.StimType]);
    Info.ds1n.isubseq   = iSubSeqs(2);
    Info.ds1n.indepval  = IndepVals(2);
    Info.ds1n.indepunit = ds1n.Stim.Presentation.X.ParUnit;
else
    Info.ds1n.iseq      = NaN;
    Info.ds1n.seqid     = NaN;
    Info.ds1n.isubseq   = NaN;
    Info.ds1n.indepval  = NaN;
    Info.ds1n.indepunit = NaN;
end

Info.ds2p.filename  = lower(ds2p.ID.Experiment.ID.Name);
Info.ds2p.icell     = ds2p.ID.iCell;
Info.ds2p.iseq      = ds2p.ID.iDataset;
Info.ds2p.seqid     = lower([num2str(ds2p.ID.iCell) '-' num2str(ds2p.ID.iRecOfCell) '-' ds2p.StimType]);
Info.ds2p.isubseq   = iSubSeqs(3);
Info.ds2p.indepval  = IndepVals(3);
Info.ds2p.indepunit = ds2p.Stim.Presentation.X.ParUnit;

if ~isvoid(ds2n)
    Info.ds2n.iseq      = ds2n.ID.iDataset;
    Info.ds2n.seqid     = lower([num2str(ds2n.ID.iCell) '-' num2str(ds2n.ID.iRecOfCell) '-' ds2n.StimType]);
    Info.ds2n.isubseq   = iSubSeqs(4);
    Info.ds2n.indepval  = IndepVals(4);
    Info.ds2n.indepunit = ds2n.Stim.Presentation.X.ParUnit;
else
    Info.ds2n.iseq      = NaN;
    Info.ds2n.seqid     = NaN;
    Info.ds2n.isubseq   = NaN;
    Info.ds2n.indepval  = NaN;
    Info.ds2n.indepunit = NaN;
end

if isnan(Info.ds1n.isubseq)
    Info.idstr1 = sprintf('%s %s#%d@%.0f%s', upper(Info.ds1p.filename), ...
        upper(Info.ds1p.seqid), Info.ds1p.iseq, Info.ds1p.indepval, Info.ds1p.indepunit);
else
    Info.idstr1 = sprintf('%s %s#%d@%.0f%s & %s#%d@%.0f%s', upper(Info.ds1p.filename), ...
        upper(Info.ds1p.seqid), Info.ds1p.iseq, Info.ds1p.indepval, Info.ds1p.indepunit, ...
        upper(Info.ds1n.seqid), Info.ds1n.iseq, Info.ds1n.indepval, Info.ds1n.indepunit);
end   
if isnan(Info.ds2n.isubseq)
    Info.idstr2 = sprintf('%s %s#%d@%.0f%s', upper(Info.ds2p.filename), ...
        upper(Info.ds2p.seqid) , Info.ds2p.iseq, Info.ds2p.indepval, Info.ds2p.indepunit);
else
    Info.idstr2 = sprintf('%s %s#%d@%.0f%s & %s#%d@%.0f%s', upper(Info.ds2p.filename), ...
        upper(Info.ds2p.seqid), Info.ds2p.iseq, Info.ds2p.indepval, Info.ds2p.indepunit, ...
        upper(Info.ds2n.seqid), Info.ds2n.iseq, Info.ds2n.indepval, Info.ds2n.indepunit);
end
Info.capstr = sprintf('%s versus %s', Info.idstr1, Info.idstr2);
Info.hdrstr = sprintf('%s \\leftrightarrow %s', Info.idstr1, Info.idstr2);

%Collecting and reorganizing stimulus parameters ...
StimParam = GetStimParam(ds1p, ds1n, ds2p, ds2n, iSubSeqs);

%Substitution of shortcuts in properties ...
if isinf(Param.anwin(2)), Param.anwin(2) = min(StimParam.burstdur); end

%Format parameter information ...
if isnan(Param.calcdf)
    CalcDFStr = 'auto';
elseif ischar(Param.calcdf)
    CalcDFStr = lower(Param.calcdf);
else
    CalcDFStr = Param2Str(Param.calcdf, 'Hz', 0);
end 
s = sprintf('AnWin = %s', Param2Str(Param.anwin, 'ms', 0));
s = char(s, sprintf('BinWidth = %s', Param2Str(Param.corbinwidth, 'ms', 2)));
s = char(s, sprintf('MaxLag = %s', Param2Str(Param.cormaxlag, 'ms', 0)));
s = char(s, sprintf('Calc. DF = %s', CalcDFStr));
s = char(s, sprintf('RunAv(Env) = %.2f(%s)', Param.envrunav, Param.envrunavunit));
s = char(s, sprintf('RunAv(Dft on COR) = %s', Param2Str(Param.corfftrunav, 'Hz', 0)));
s = char(s, sprintf('RunAv(Dft on DIF) = %s', Param2Str(Param.diffftrunav, 'Hz', 0)));
Param.str = s;

%----------------------------------------------------------------------------
function CheckParam(Param)

%Syntax parameters ...
if ~any(strcmpi(Param.subseqinput, {'indepval', 'subseq'})), error('Property subseqinput must be ''indepval'' or ''subseq''.'); end

%Calculation parameters ...
if ~isnumeric(Param.anwin) | (size(Param.anwin) ~= [1,2]) | ~isinrange(Param.anwin, [0, +Inf]), error('Invalid value for property anwin.'); end
if ~isnumeric(Param.corbinwidth) || (length(Param.corbinwidth) ~= 1) || (Param.corbinwidth <= 0), error('Invalid value for property corbinwidth.'); end
if ~isnumeric(Param.cormaxlag) || (length(Param.cormaxlag) ~= 1) || (Param.cormaxlag <= 0), error('Invalid value for property cormaxlag.'); end
if ~any(strcmpi(Param.envrunavunit, {'#', 'ms'})), error('Property envrunavunit must be ''#'' or ''ms''.'); end
if ~isnumeric(Param.envrunav) || (length(Param.envrunav) ~= 1) || (Param.envrunav < 0), error('Invalid value for property envrunav.'); end
if ~isnumeric(Param.corfftrunav) || (length(Param.corfftrunav) ~= 1) || (Param.corfftrunav < 0), error('Invalid value for property corfftrunav.'); end
if ~isnumeric(Param.diffftrunav) || (length(Param.diffftrunav) ~= 1) || (Param.diffftrunav < 0), error('Invalid value for property diffftrunav.'); end
if ~isnumeric(Param.sumfftrunav) || (length(Param.sumfftrunav) ~= 1) || (Param.sumfftrunav < 0), error('Invalid value for property diffftrunav.'); end
if ~(isnumeric(Param.calcdf) && ((Param.calcdf > 0) || isnan(Param.calcdf))) && ...
        ~(ischar(Param.calcdf) && any(strcmpi(Param.calcdf, {'cf', 'df'})))
    error('Property calcdf must be positive integer, NaN, ''cf'' or ''df''.'); 
end
if ~any(strncmpi(Param.calctype, {'scc', 'avgsac'}, 1)) 
    error('Property calctype must be ''scc'' or ''avgsac''.'); 
end

%Plot parameters ...
if ~any(strcmpi(Param.plot, {'yes', 'no'})), error('Property plot must be ''yes'' or ''no''.'); end
if ~isinrange(Param.corxrange, [-Inf +Inf]), error('Invalid value for property corxrange.'); end
if ~isnumeric(Param.corxstep) || (length(Param.corxstep) ~= 1) || (Param.corxstep <= 0), error('Invalid value for property corxstep.'); end
if ~isinrange(Param.fftxrange, [0 +Inf]), error('Invalid value for property fftxrange.'); end
if ~isnumeric(Param.fftxstep) || (length(Param.fftxstep) ~= 1) || (Param.fftxstep <= 0), error('Invalid value for property fftxstep.'); end
if ~any(strcmpi(Param.fftyunit, {'dB', 'P'})), error('Property fftyunit must be ''dB'' or ''P''.'); end
if ~isinrange(Param.fftyrange, [-Inf +Inf]), error('Invalid value for property fftyrange.'); end

%----------------------------------------------------------------------------
function StimParam = GetStimParam(ds1p, ds1n, ds2p, ds2n, iSubSeqs)

dsNames = {'ds1p', 'ds1n', 'ds2p', 'ds2n'}; Nds = 4;

[StimParam.burstdur, StimParam.repdur, StimParam.nrep] = deal(NaN(1, Nds));
StimParam.spl = NaN(2, Nds);
for n = 1:Nds
    ds = eval(dsNames{n});
    if ~isnan(iSubSeqs(n))
        %If stimulus or repetition duration are not the same for different channels then
        %the minimum value is used ...
        StimParam.burstdur(n) = round(min(ds.Stim.BurstDur));
        StimParam.repdur(n) = round(min(ds.Stim.ISI));
        StimParam.nrep(n) = round(min(ds.Stim.Presentation.Nrep));
        SPL = GetSPL(ds); StimParam.spl(:, n) = round(SPL(iSubSeqs(n), [1, end])');
    end
end
StimParam.avgspl1 = CombineSPLs(denan(StimParam.spl(:, [1 2]))');
StimParam.avgspl2 = CombineSPLs(denan(StimParam.spl(:, [3 4]))');
StimParam.avgspl  = CombineSPLs(denan(StimParam.spl(:))');

%Format stimulus parameters ...
s = sprintf('BurstDur = %s ms', mat2str(StimParam.burstdur));
s = char(s, sprintf('IntDur = %s ms', mat2str(StimParam.repdur)));
s = char(s, sprintf('#Reps = %s', mat2str(StimParam.nrep)));
s = char(s, sprintf('SPL = %s dB', mat2str(CombineSPLs(StimParam.spl')')));
StimParam.str = s;

%----------------------------------------------------------------------------
function Str = Param2Str(V, Unit, Prec)

C = num2cell(V);
Sz = size(V);
N  = prod(Sz);

if (N == 1) || all(isequal(C{:}))
    Str = sprintf(['%.'  int2str(Prec) 'f%s'], V(1), Unit);
elseif (N == 2)
    Str = sprintf(['%.' int2str(Prec) 'f%s/%.' int2str(Prec) 'f%s'], V(1), Unit, V(2), Unit);
elseif any(Sz == 1)
    Str = sprintf(['%.' int2str(Prec) 'f%s..%.' int2str(Prec) 'f%s'], min(V(:)), Unit, max(V(:)), Unit); 
else
    Str = sprintf(['%.' int2str(Prec) 'f%s/%.' int2str(Prec) 'f%s..%.' int2str(Prec) 'f%s/%.' int2str(Prec) 'f%s'], ...
        min(V(:, 1)), Unit, min(V(:, 2)), Unit, max(V(:, 1)), Unit, max(V(:, 2)), Unit); 
end

%----------------------------------------------------------------------------
function RC = CalcRC(SptP, SptN, Param)

WinDur = abs(diff(Param.anwin));

Rp = 1e3*mean(cellfun('length', SptP))/WinDur;
if ~isempty(SptN), Rn = 1e3*mean(cellfun('length', SptN))/WinDur; else, Rn = []; end

RC.mean = mean([Rp, Rn]);
RC.str = sprintf('AvgR = %s', Param2Str(RC.mean, 'spk/sec', 0));

%----------------------------------------------------------------------------
function [SXAC, DAC, MAC] = CalcAC(SptP, SptN, Thr, Param)

WinDur = abs(diff(Param.anwin)); %Duration of analysis window in ms ...
SptP = anwin(SptP, Param.anwin);
SptN = anwin(SptN, Param.anwin);

if ~isempty(SptN)
    %Correlation of noise token A+ responses of a cell with the responses of that same cell to that same noise
    %token. If spiketrains are derived from the same cell this is called a Shuffled Auto-
    %Correlogram (or SAC). 'Shuffled' because of the shuffling of repetitions in order to avoid to correlation
    %of a repetition with itself. The terminolgy AutoCorrelogram is only used when comparing spiketrains 
    %collected from the same cell.
    [Ypp, ~, NC] = SPTCORR(SptP, 'nodiag', Param.cormaxlag, Param.corbinwidth, WinDur); %SAC ...
    Ypp = ApplyNorm(Ypp, NC);
    %Correlation of noise token A- responses of a cell with the responses of that same cell to that same noise
    %token.
    [Ynn, ~, NC] = SPTCORR(SptN, 'nodiag', Param.cormaxlag, Param.corbinwidth, WinDur); %SAC ...
    Ynn = ApplyNorm(Ynn, NC);
    %Correlation of noise token A+ responses of a cell with the responses of that same cell to a different noise
    %token, in this case A-. Because of the fact that we correlate across stimuli this type of correlogram is 
    %designated XAC, when comparing responses from the same cell.
    [Ypn, ~, NC] = SPTCORR(SptP, SptN, Param.cormaxlag, Param.corbinwidth, WinDur); %XAC ...
    Ypn = ApplyNorm(Ypn, NC);
    %Correlation of noise token A- responses of a cell with the responses of that same cell to a different noise
    %token, in this case A+.
    [Ynp, T, NC] = SPTCORR(SptN, SptP, Param.cormaxlag, Param.corbinwidth, WinDur); %XAC ...
    Ynp = ApplyNorm(Ynp, NC);
    
    %Calculation of the DIFCOR by taking the average of the two SACs and the two XACs and subtracting the second
    %from the first ...
    Ysac = mean([Ypp; Ynn]); 
    Yxac = mean([Ypn; Ynp]);
    Ydifcor = Ysac - Yxac;
    Ysumcor = Ysac + Yxac;
    
    %Performing spectrum analysis on the DIFCOR and SUMCOR. Because a difcor has no DC component in comparison with
    %other correlograms, this almost always results in a representative dominant frequency ...
    FFTdif = spectana(T, Ydifcor, 'RunAvUnit', 'Hz', 'RunAvRange', Param.diffftrunav);
    FFTsum = spectana(T, Ysumcor, 'RunAvUnit', 'Hz', 'RunAvRange', Param.sumfftrunav);
    %The magnitude spectrum of a correlogram function is actually a power spectrum, therefore all
    %magnitude units need to be changed ...
    FFTdif.Magn.P  = FFTdif.Magn.A;
    FFTdif.Magn.A  = sqrt(FFTdif.Magn.A);
    FFTdif.Magn.dB = FFTdif.Magn.dB/2;
    FFTsum.Magn.P  = FFTsum.Magn.A;
    FFTsum.Magn.A  = sqrt(FFTsum.Magn.A);
    FFTsum.Magn.dB = FFTsum.Magn.dB/2;
    
    %Performing spectrum analysis on the SAC. Because an autocorrelogram has a DC component this is
    %removed first ...
    FFTsac = spectana(T, detrend(Ysac, 'constant'), 'RunAvUnit', 'Hz', 'RunAvRange', Param.corfftrunav);
    %The magnitude spectrum of a correlogram function is actually a power spectrum, therefore all
    %magnitude units need to be changed ...
    FFTsac.Magn.P  = FFTsac.Magn.A;
    FFTsac.Magn.A  = sqrt(FFTsac.Magn.A);
    FFTsac.Magn.dB = FFTsac.Magn.dB/2;
    
    %Determine which dominant frequency to be used in the calculation ...
    DomFreq = DetermineCalcDF(Param.calcdf, Thr.cf, FFTdif.DF, FFTsac.DF);
    %Dominant period in ms ...
    if (DomFreq ~= 0)
        DomPer = 1000/DomFreq;
    else
        DomPer = NaN;
    end
    
    %Calculating the half height width on the peak of the SAC (relative to asymptote one)...
    HalfMaxSac = ((max(Ysac)-1)/2)+1;
    SacHHWx = cintersect(T, Ysac, HalfMaxSac); SacHHW = abs(diff(SacHHWx));
    
    %Calculating envelope and Half Height Width of DIFCOR and SUMCOR ...
    if strcmpi(Param.envrunavunit, 'ms')
        EnvRunAvN = round(Param.envrunav/Param.corbinwidth);
    else
        EnvRunAvN = round((Param.envrunav*DomPer)/Param.corbinwidth);
    end
    YenvDif = runav(abs(hilbert(Ydifcor)), EnvRunAvN); 
    HalfMaxEnvDif = max(YenvDif)/2;
    DifHHWx = cintersect(T, YenvDif, HalfMaxEnvDif); DifHHW = abs(diff(DifHHWx));
    YenvSum = runav(abs(hilbert(Ysumcor)), EnvRunAvN); 
    HalfMaxEnvSum = max(YenvSum)/2;
    SumHHWx = cintersect(T, YenvSum, HalfMaxEnvSum); SumHHW = abs(diff(SumHHWx));
    
    %Extracting information on secundary peaks in the DIFCOR and its enveloppe ...
    [~, ~, DifXsecPeaks, DifYsecPeaks] = getPeaks(T, Ydifcor, 0, DomPer);
    [~, ~, SumXsecPeaks, SumYsecPeaks] = getPeaks(T, Ysumcor, 0, DomPer);
    
    %Reorganizing calculated data ...
    SXAC.lag      = T;
    SXAC.normco   = [Ysac; Yxac];
    SXAC.max      = max(Ysac);
    SXAC.hhw      = SacHHW;
    SXAC.hhwx     = SacHHWx;
    SXAC.halfmax  = HalfMaxSac;
    SXAC.fft.freq = FFTsac.Freq;
    SXAC.fft.p    = FFTsac.Magn.P;   
    SXAC.fft.db   = FFTsac.Magn.dB;
    SXAC.fft.df   = FFTsac.DF;
    SXAC.fft.bw   = FFTsac.BW;
    
    DAC.lag         = T;
    DAC.normco      = [Ydifcor; YenvDif; -YenvDif];
    DAC.max         = max(Ydifcor);
    DAC.maxsecpks   = DifYsecPeaks;
    DAC.lagsecpks   = DifXsecPeaks;
    DAC.fft.freq    = FFTdif.Freq;
    DAC.fft.p       = FFTdif.Magn.P;   
    DAC.fft.db      = FFTdif.Magn.dB;
    DAC.fft.df      = FFTdif.DF;
    DAC.fft.bw      = FFTdif.BW;
    DAC.env.hhw     = DifHHW;
    DAC.env.hhwx    = DifHHWx;
    DAC.env.halfmax = HalfMaxEnvDif;
    
    MAC.lag         = T;
    MAC.normco      = [Ysumcor; YenvSum; -YenvSum];
    MAC.max         = max(Ysumcor);
    MAC.maxsecpks   = SumYsecPeaks;
    MAC.lagsecpks   = SumXsecPeaks;
    MAC.fft.freq    = FFTsum.Freq;
    MAC.fft.p       = FFTsum.Magn.P;   
    MAC.fft.db      = FFTsum.Magn.dB;
    MAC.fft.df      = FFTsum.DF;
    MAC.fft.bw      = FFTsum.BW;
    MAC.env.hhw     = SumHHW;
    MAC.env.hhwx    = SumHHWx;
    MAC.env.halfmax = HalfMaxEnvSum;
else
    %Correlation of noise token A+ responses of a cell with the responses of that same cell to that same noise
    %token ...
    [Ysac, T, NC] = SPTCORR(SptP, 'nodiag', Param.cormaxlag, Param.corbinwidth, WinDur); %SAC ...
    Ysac = ApplyNorm(Ysac, NC);
    
    %Performing spectrum analysis on the SAC. Because an autocorrelogram has a DC component this is
    %removed first ...
    FFTsac = spectana(T, detrend(Ysac, 'constant'), 'RunAvUnit', 'Hz', 'RunAvRange', Param.corfftrunav);
    %The magnitude spectrum of a correlogram function is actually a power spectrum, therefore all
    %magnitude units need to be changed ...
    FFTsac.Magn.P  = FFTsac.Magn.A;
    FFTsac.Magn.A  = sqrt(FFTsac.Magn.A);
    FFTsac.Magn.dB = FFTsac.Magn.dB/2;
    
    %Determine which dominant frequency to be used in the calculation ...
    DomFreq = DetermineCalcDF(Param.calcdf, Thr.cf, NaN, FFTsac.DF);
    %Dominant period in ms ...
    if (DomFreq ~= 0)
        DomPer = 1000/DomFreq;
    else
        DomPer = NaN;
    end
    
    %Calculating the half height width on the peak of the SAC (relative to asymptote one)...
    HalfMaxSac = ((max(Ysac)-1)/2)+1;
    SacHHWx = cintersect(T, Ysac, HalfMaxSac); SacHHW = abs(diff(SacHHWx));
    
    %Reorganizing calculated data ...
    SXAC.lag      = T;
    SXAC.normco   = Ysac;
    SXAC.max      = max(Ysac);
    SXAC.hhw      = SacHHW;
    SXAC.hhwx     = SacHHWx;
    SXAC.halfmax  = HalfMaxSac;
    SXAC.fft.freq = FFTsac.Freq;
    SXAC.fft.p    = FFTsac.Magn.P;   
    SXAC.fft.db   = FFTsac.Magn.dB;
    SXAC.fft.df   = FFTsac.DF;
    SXAC.fft.bw   = FFTsac.BW;
    
    DAC.lag         = [];
    DAC.normco      = [];
    DAC.max         = NaN;
    DAC.maxsecpks   = [NaN, NaN];
    DAC.lagsecpks   = [NaN, NaN];
    DAC.fft.freq    = [];
    DAC.fft.p       = [];   
    DAC.fft.db      = [];
    DAC.fft.df      = NaN;
    DAC.fft.bw      = NaN;
    DAC.env.hhw     = NaN;
    DAC.env.hhwx    = [NaN, NaN];
    DAC.env.halfmax = NaN;
end

%----------------------------------------------------------------------------
function [SXCC, DCC, MCC] = CalcAvgAC(Spt1P, Spt1N, Spt2P, Spt2N, Thr, Param)

WinDur = abs(diff(Param.anwin)); %Duration of analysis window in ms ...
Spt1P = ANWIN(Spt1P, Param.anwin);
Spt1N = ANWIN(Spt1N, Param.anwin);
Spt2P = ANWIN(Spt2P, Param.anwin);
Spt2N = ANWIN(Spt2N, Param.anwin);

if ~isempty(Spt1N) && ~isempty(Spt2N)
    [Ypp1NCo, ~, NC] = SPTCORR(Spt1P, 'nodiag', Param.cormaxlag, Param.corbinwidth, WinDur); %SAC ...
    Ypp1 = ApplyNorm(Ypp1NCo, NC); Ypp1Rate = ApplyNorm(Ypp1NCo, NC, 'rate');

    [Ynn1NCo, ~, NC] = SPTCORR(Spt1N, 'nodiag', Param.cormaxlag, Param.corbinwidth, WinDur); %SAC ...
    Ynn1 = ApplyNorm(Ynn1NCo, NC); Ynn1Rate = ApplyNorm(Ynn1NCo, NC, 'rate');
    
    [Ypp2NCo, ~, NC] = SPTCORR(Spt2P, 'nodiag', Param.cormaxlag, Param.corbinwidth, WinDur); %SAC ...
    Ypp2 = ApplyNorm(Ypp2NCo, NC); Ypp2Rate = ApplyNorm(Ypp2NCo, NC, 'rate');

    [Ynn2NCo, ~, NC] = SPTCORR(Spt2N, 'nodiag', Param.cormaxlag, Param.corbinwidth, WinDur); %SAC ...
    Ynn2 = ApplyNorm(Ynn2NCo, NC); Ynn2Rate = ApplyNorm(Ynn2NCo, NC, 'rate');
    
    [Ypn1NCo, ~, NC] = SPTCORR(Spt1P, Spt1N, Param.cormaxlag, Param.corbinwidth, WinDur); %XAC ...
    Ypn1 = ApplyNorm(Ypn1NCo, NC); Ypn1Rate = ApplyNorm(Ypn1NCo, NC, 'rate');

    [Ynp1NCo, ~, NC] = SPTCORR(Spt1N, Spt1P, Param.cormaxlag, Param.corbinwidth, WinDur); %XAC ...
    Ynp1 = ApplyNorm(Ynp1NCo, NC); Ynp1Rate = ApplyNorm(Ynp1NCo, NC, 'rate');

    [Ypn2NCo, ~, NC] = SPTCORR(Spt2P, Spt2N, Param.cormaxlag, Param.corbinwidth, WinDur); %XAC ...
    Ypn2 = ApplyNorm(Ypn2NCo, NC); Ypn2Rate = ApplyNorm(Ypn2NCo, NC, 'rate');

    [Ynp2NCo, T, NC] = SPTCORR(Spt2N, Spt2P, Param.cormaxlag, Param.corbinwidth, WinDur); %XAC ...
    Ynp2 = ApplyNorm(Ynp2NCo, NC); Ynp2Rate = ApplyNorm(Ynp2NCo, NC, 'rate');

    Ysac = mean([Ypp1; Ynn1; Ypp2; Ynn2]); YsacRate = mean([Ypp1Rate; Ynn1Rate; Ypp2Rate; Ynn2Rate]);
    Yxac = mean([Ypn1; Ynp1; Ypn2; Ynp2]); YxacRate = mean([Ypn1Rate; Ynp1Rate; Ypn2Rate; Ynp2Rate]);
    Ydifcor = Ysac - Yxac; YdifcorRate = YsacRate - YxacRate;
    Ysumcor = Ysac + Yxac; YsumcorRate = YsacRate + YxacRate;
    
    %Performing spectrum analysis on the DIFCOR and SUMCOR. Because a difcor has no DC component in comparison with
    %other correlograms, this almost always results in a representative dominant frequency ...
    FFTdif = spectana(T, Ydifcor, 'RunAvUnit', 'Hz', 'RunAvRange', Param.diffftrunav);
    FFTsum = spectana(T, Ysumcor, 'RunAvUnit', 'Hz', 'RunAvRange', Param.sumfftrunav);
    %The magnitude spectrum of a correlogram function is actually a power spectrum, therefore all
    %magnitude units need to be changed ...
    FFTdif.Magn.P  = FFTdif.Magn.A;
    FFTdif.Magn.A  = sqrt(FFTdif.Magn.A);
    FFTdif.Magn.dB = FFTdif.Magn.dB/2;
    FFTsum.Magn.P  = FFTsum.Magn.A;
    FFTsum.Magn.A  = sqrt(FFTsum.Magn.A);
    FFTsum.Magn.dB = FFTsum.Magn.dB/2;
    
    %Performing spectrum analysis on the SAC. Because a crosscorrelogram has a DC component this is
    %removed first ...
    FFTsac = spectana(T, detrend(Ysac, 'constant'), 'RunAvUnit', 'Hz', 'RunAvRange', Param.corfftrunav);
    %The magnitude spectrum of a correlogram function is actually a power spectrum, therefore all
    %magnitude units need to be changed ...
    FFTsac.Magn.P  = FFTsac.Magn.A;
    FFTsac.Magn.A  = sqrt(FFTsac.Magn.A);
    FFTsac.Magn.dB = FFTsac.Magn.dB/2;
    
    %Determine which dominant frequency to be used in the calculation ...
    DomFreq = DetermineCalcDF(Param.calcdf, mean(deNaN(cat(2, Thr.cf))), FFTdif.DF, FFTsac.DF);
    %Dominant period in ms ...
    if (DomFreq ~= 0)
        DomPer = 1000/DomFreq;
    else
        DomPer = NaN;
    end
    
    %Calculating the half height width on the peak of the SAC (relative to asymptote one)...
    HalfMaxSac = ((max(Ysac)-1)/2)+1;
    SacHHWx = cintersect(T, Ysac, HalfMaxSac); SacHHW = abs(diff(SacHHWx));
    
    %Extracting information on secundary peaks in the SCC ...
    [~, ~, SacXsecPeaks, SacYsecPeaks] = getPeaks(T, Ysac, 0, DomPer);
    
    %Calculating envelope and Half Height Width of DIFCOR and SUMCOR ...
    if strcmpi(Param.envrunavunit, 'ms')
        EnvRunAvN = round(Param.envrunav/Param.corbinwidth);
    else
        EnvRunAvN = round((Param.envrunav*DomPer)/Param.corbinwidth);
    end
    YenvDif = runav(abs(hilbert(Ydifcor)), EnvRunAvN); 
    HalfMaxEnvDif = max(YenvDif)/2;
    DifHHWx = cintersect(T, YenvDif, HalfMaxEnvDif); DifHHW = abs(diff(DifHHWx));
    YenvSum = runav(abs(hilbert(Ysumcor)), EnvRunAvN); 
    HalfMaxEnvSum = max(YenvSum)/2;
    SumHHWx = cintersect(T, YenvSum, HalfMaxEnvSum); SumHHW = abs(diff(SumHHWx));
    
    %Extracting information on secundary peaks in the DIFCOR and SUMCOR ...
    [~, ~, DifXsecPeaks, DifYsecPeaks] = getPeaks(T, Ydifcor, 0, DomPer);
    [~, ~, SumXsecPeaks, SumYsecPeaks] = getPeaks(T, Ysumcor, 0, DomPer);
    
    %Reorganizing calculated data ...
    SXCC.lag        = T;
    SXCC.normco     = [Ysac; Yxac];
    [SXCC.max, idx] = max(Ysac);
    SXCC.rate       = max(YsacRate);
    SXCC.lagatmax   = T(idx);
    SXCC.maxsecpks  = SacYsecPeaks;
    SXCC.lagsecpks  = SacXsecPeaks;
    SXCC.hhw        = SacHHW;
    SXCC.hhwx       = SacHHWx;
    SXCC.halfmax    = HalfMaxSac;
    SXCC.fft.freq   = FFTsac.Freq;
    SXCC.fft.p      = FFTsac.Magn.P;   
    SXCC.fft.db     = FFTsac.Magn.dB;
    SXCC.fft.df     = FFTsac.DF;
    SXCC.fft.bw     = FFTsac.BW;
    
    DCC.lag          = T;
    DCC.normco       = [Ydifcor; YenvDif; -YenvDif];
    [DCC.max, idx]   = max(Ydifcor);
    DCC.rate         = max(YdifcorRate);
    DCC.lagatmax     = T(idx);
    DCC.maxsecpks    = DifYsecPeaks;
    DCC.lagsecpks    = DifXsecPeaks;
    DCC.fft.freq     = FFTdif.Freq;
    DCC.fft.p        = FFTdif.Magn.P;   
    DCC.fft.db       = FFTdif.Magn.dB;
    DCC.fft.df       = FFTdif.DF;
    DCC.fft.bw       = FFTdif.BW;
    [DCC.env.max, idx]= max(YenvDif);
    DCC.env.lagatmax = T(idx);
    DCC.env.hhw      = DifHHW;
    DCC.env.hhwx     = DifHHWx;
    DCC.env.halfmax  = HalfMaxEnvDif;
    
    MCC.lag          = T;
    MCC.normco       = [Ysumcor; YenvSum; -YenvSum];
    [MCC.max, idx]   = max(Ysumcor);
    MCC.rate         = max(YsumcorRate);
    MCC.lagatmax     = T(idx);
    MCC.maxsecpks    = SumYsecPeaks;
    MCC.lagsecpks    = SumXsecPeaks;
    MCC.fft.freq     = FFTsum.Freq;
    MCC.fft.p        = FFTsum.Magn.P;   
    MCC.fft.db       = FFTsum.Magn.dB;
    MCC.fft.df       = FFTsum.DF;
    MCC.fft.bw       = FFTsum.BW;
    [MCC.env.max, idx]= max(YenvSum);
    MCC.env.lagatmax = T(idx);
    MCC.env.hhw      = SumHHW;
    MCC.env.hhwx     = SumHHWx;
    MCC.env.halfmax  = HalfMaxEnvSum;
else
    warning(['Calculation of average SAC ony possible if responses of' ...
        'second fiber are supplied,\n or more appropriate for this' ...
        'analysis the response of the same fiber to different noise token.']);

    SXCC.lag        = [];
    SXCC.normco     = [];
    SXCC.max        = NaN;
    SXCC.rate       = NaN;
    SXCC.lagatmax   = NaN;
    SXCC.maxsecpks  = [NaN, NaN];
    SXCC.lagsecpks  = [NaN, NaN];
    SXCC.hhw        = NaN;
    SXCC.hhwx       = [NaN, NaN];
    SXCC.halfmax    = NaN;
    SXCC.fft.freq   = NaN;
    SXCC.fft.p      = NaN;   
    SXCC.fft.db     = NaN;
    SXCC.fft.df     = NaN;
    SXCC.fft.bw     = NaN;
    
    DCC.lag          = [];
    DCC.normco       = [];
    DCC.max          = NaN;
    DCC.rate         = NaN;
    DCC.lagatmax     = NaN;
    DCC.maxsecpks    = [NaN, NaN];
    DCC.lagsecpks    = [NaN, NaN];
    DCC.fft.freq     = [];
    DCC.fft.p        = [];   
    DCC.fft.db       = [];
    DCC.fft.df       = NaN;
    DCC.fft.bw       = NaN;
    DCC.env.max      = NaN;
    DCC.env.lagatmax = NaN;
    DCC.env.hhw      = NaN;
    DCC.env.hhwx     = [NaN, NaN];
    DCC.env.halfmax  = NaN;
end

%----------------------------------------------------------------------------
function [SXCC, DCC, MCC] = CalcCC(Spt1P, Spt1N, Spt2P, Spt2N, Thr, Param)

WinDur = abs(diff(Param.anwin)); %Duration of analysis window in ms ...
Spt1P = anwin(Spt1P, Param.anwin);
Spt1N = anwin(Spt1N, Param.anwin);
Spt2P = anwin(Spt2P, Param.anwin);
Spt2N = anwin(Spt2N, Param.anwin);

if ~isempty(Spt1N) && ~isempty(Spt2N)
    %Correlation of noise token A+ responses of a cell with the responses of another cell to that same noise
    %token. If spiketrains are derived from different cells this is called a Shuffled Cross-
    %Correlogram (or SCC). 'Shuffled' because of the shuffling of repetitions in order to avoid to correlation
    %of a repetition with itself. The terminolgy CrossCorrelogram is only used when comparing spiketrains 
    %collected from different cells.
    [YppNCo, ~, NC] = SPTCORR(Spt1P, Spt2P, Param.cormaxlag, Param.corbinwidth, WinDur); %SCC ...
    Ypp = ApplyNorm(YppNCo, NC);
    YppRate = ApplyNorm(YppNCo, NC, 'rate');
    %Correlation of noise token A- responses of a cell with the responses of a different cell to that same noise
    %token.
    [YnnNCo, ~, NC] = SPTCORR(Spt1N, Spt2N, Param.cormaxlag, Param.corbinwidth, WinDur); %SCC ...
    Ynn = ApplyNorm(YnnNCo, NC);
    YnnRate = ApplyNorm(YnnNCo, NC, 'rate');
    %Correlation of noise token A+ responses of a cell with the responses of a different cell to a different noise
    %token, in this case A-. Because of the fact that we correlate across stimuli this type of correlogram is 
    %designated XCC.
    [YpnNCo, ~, NC] = SPTCORR(Spt1P, Spt2N, Param.cormaxlag, Param.corbinwidth, WinDur); %XCC ...
    Ypn = ApplyNorm(YpnNCo, NC);
    YpnRate = ApplyNorm(YpnNCo, NC, 'rate');
    %Correlation of noise token A- responses of a cell with the responses of a different cell to a different noise
    %token, in this case A+.
    [YnpNCo, T, NC] = SPTCORR(Spt1N, Spt2P, Param.cormaxlag, Param.corbinwidth, WinDur); %XCC ...
    Ynp = ApplyNorm(YnpNCo, NC);
    YnpRate = ApplyNorm(YnpNCo, NC, 'rate');
    
    %Calculation of the DIFCOR by taking the average of the two SCCs and the two XCCs and subtracting the second
    %from the first ...
    Yscc = mean([Ypp; Ynn]);
    YsccRate = mean([YppRate; YnnRate]);
    Yxcc = mean([Ypn; Ynp]);
    YxccRate = mean([YpnRate; YnpRate]);
    Ydifcor = Yscc - Yxcc;
    YdifcorRate = YsccRate - YxccRate;
    Ysumcor = Yscc + Yxcc;
    YsumcorRate = YsccRate + YxccRate;
    
    %Performing spectrum analysis on the DIFCOR and SUMCOR. Because a difcor has no DC component in comparison with
    %other correlograms, this almost always results in a representative dominant frequency ...
    FFTdif = spectana(T, Ydifcor, 'RunAvUnit', 'Hz', 'RunAvRange', Param.diffftrunav);
    FFTsum = spectana(T, Ysumcor, 'RunAvUnit', 'Hz', 'RunAvRange', Param.sumfftrunav);
    %The magnitude spectrum of a correlogram function is actually a power spectrum, therefore all
    %magnitude units need to be changed ...
    FFTdif.Magn.P  = FFTdif.Magn.A;
    FFTdif.Magn.A  = sqrt(FFTdif.Magn.A);
    FFTdif.Magn.dB = FFTdif.Magn.dB/2;
    FFTsum.Magn.P  = FFTsum.Magn.A;
    FFTsum.Magn.A  = sqrt(FFTsum.Magn.A);
    FFTsum.Magn.dB = FFTsum.Magn.dB/2;
    
    %Performing spectrum analysis on the SCC. Because a crosscorrelogram has a DC component this is
    %removed first ...
    FFTscc = spectana(T, detrend(Yscc, 'constant'), 'RunAvUnit', 'Hz', 'RunAvRange', Param.corfftrunav);
    %The magnitude spectrum of a correlogram function is actually a power spectrum, therefore all
    %magnitude units need to be changed ...
    FFTscc.Magn.P  = FFTscc.Magn.A;
    FFTscc.Magn.A  = sqrt(FFTscc.Magn.A);
    FFTscc.Magn.dB = FFTscc.Magn.dB/2;
    
    %Determine which dominant frequency to be used in the calculation ...
    DomFreq = DetermineCalcDF(Param.calcdf, mean(denan(cat(2, Thr.cf))), FFTdif.DF, FFTscc.DF);
    %Dominant period in ms ...
    if (DomFreq ~= 0)
        DomPer = 1000/DomFreq;
    else
        DomPer = NaN;
    end
    
    %Calculating the half height width on the peak of the SCC (relative to asymptote one)...
    HalfMaxScc = ((max(Yscc)-1)/2)+1;
    SccHHWx = cintersect(T, Yscc, HalfMaxScc); SccHHW = abs(diff(SccHHWx));
    
    %Extracting information on secundary peaks in the SCC ...
    [~, ~, SccXsecPeaks, SccYsecPeaks] = getPeaks(T, Yscc, 0, DomPer);
    
    %Calculating envelope and Half Height Width of DIFCOR and SUMCOR ...
    if strcmpi(Param.envrunavunit, 'ms')
        EnvRunAvN = round(Param.envrunav/Param.corbinwidth);
    else
        EnvRunAvN = round((Param.envrunav*DomPer)/Param.corbinwidth);
    end
    YenvDif = runav(abs(hilbert(Ydifcor)), EnvRunAvN); 
    HalfMaxEnvDif = max(YenvDif)/2;
    DifHHWx = cintersect(T, YenvDif, HalfMaxEnvDif); DifHHW = abs(diff(DifHHWx));
    YenvSum = runav(abs(hilbert(Ysumcor)), EnvRunAvN); 
    HalfMaxEnvSum = max(YenvSum)/2;
    SumHHWx = cintersect(T, YenvSum, HalfMaxEnvSum); SumHHW = abs(diff(SumHHWx));
    
    %Extracting information on secundary peaks in the DIFCOR and SUMCOR ...
    [~, ~, DifXsecPeaks, DifYsecPeaks] = getPeaks(T, Ydifcor, 0, DomPer);
    [~, ~, SumXsecPeaks, SumYsecPeaks] = getPeaks(T, Ysumcor, 0, DomPer);
    
    %Reorganizing calculated data ...
    SXCC.lag        = T;
    SXCC.normco     = [Yscc; Yxcc];
    [SXCC.max, idx] = max(Yscc);
    SXCC.rate       = max(YsccRate);
    SXCC.lagatmax   = T(idx);
    SXCC.maxsecpks  = SccYsecPeaks;
    SXCC.lagsecpks  = SccXsecPeaks;
    SXCC.hhw        = SccHHW;
    SXCC.hhwx       = SccHHWx;
    SXCC.halfmax    = HalfMaxScc;
    SXCC.fft.freq   = FFTscc.Freq;
    SXCC.fft.p      = FFTscc.Magn.P;   
    SXCC.fft.db     = FFTscc.Magn.dB;
    SXCC.fft.df     = FFTscc.DF;
    SXCC.fft.bw     = FFTscc.BW;
    
    DCC.lag          = T;
    DCC.normco       = [Ydifcor; YenvDif; -YenvDif];
    [DCC.max, idx]   = max(Ydifcor);
    DCC.rate         = max(YdifcorRate);
    DCC.lagatmax     = T(idx);
    DCC.maxsecpks    = DifYsecPeaks;
    DCC.lagsecpks    = DifXsecPeaks;
    DCC.fft.freq     = FFTdif.Freq;
    DCC.fft.p        = FFTdif.Magn.P;   
    DCC.fft.db       = FFTdif.Magn.dB;
    DCC.fft.df       = FFTdif.DF;
    DCC.fft.bw       = FFTdif.BW;
    [DCC.env.max, idx]= max(YenvDif);
    DCC.env.lagatmax = T(idx);
    DCC.env.hhw      = DifHHW;
    DCC.env.hhwx     = DifHHWx;
    DCC.env.halfmax  = HalfMaxEnvDif;
    
    MCC.lag          = T;
    MCC.normco       = [Ysumcor; YenvSum; -YenvSum];
    [MCC.max, idx]   = max(Ysumcor);
    MCC.rate         = max(YsumcorRate);
    MCC.lagatmax     = T(idx);
    MCC.maxsecpks    = SumYsecPeaks;
    MCC.lagsecpks    = SumXsecPeaks;
    MCC.fft.freq     = FFTsum.Freq;
    MCC.fft.p        = FFTsum.Magn.P;   
    MCC.fft.db       = FFTsum.Magn.dB;
    MCC.fft.df       = FFTsum.DF;
    MCC.fft.bw       = FFTsum.BW;
    [MCC.env.max, idx]= max(YenvSum);
    MCC.env.lagatmax = T(idx);
    MCC.env.hhw      = SumHHW;
    MCC.env.hhwx     = SumHHWx;
    MCC.env.halfmax  = HalfMaxEnvSum;
else
    %Correlation of noise token A+ responses of a cell with the responses of another cell to that same noise
    %token.
    [YsccNCo, T, NC] = SPTCORR(Spt1P, Spt2P, Param.cormaxlag, Param.corbinwidth, WinDur); %SCC ...
    Yscc = ApplyNorm(YsccNCo, NC); YsccRate = ApplyNorm(YsccNCo, NC, 'Rate');

    %Performing spectrum analysis on the SCC. Because a crosscorrelogram has a DC component this is
    %removed first ...
    FFTscc = spectana(T, detrend(Yscc, 'constant'), 'RunAvUnit', 'Hz', 'RunAvRange', Param.corfftrunav);
    %The magnitude spectrum of a correlogram function is actually a power spectrum, therefore all
    %magnitude units need to be changed ...
    FFTscc.Magn.P  = FFTscc.Magn.A;
    FFTscc.Magn.A  = sqrt(FFTscc.Magn.A);
    FFTscc.Magn.dB = FFTscc.Magn.dB/2;
    
    %Determine which dominant frequency to be used in the calculation ...
    DomFreq = DetermineCalcDF(Param.calcdf, mean(denan(cat(2, Thr.cf))), NaN, FFTscc.DF);
    %Dominant period in ms ...
    if (DomFreq ~= 0)
        DomPer = 1000/DomFreq;
    else
        DomPer = NaN;
    end
    
    %Calculating the half height width on the peak of the SCC (relative to asymptote one)...
    HalfMaxScc = ((max(Yscc)-1)/2)+1;
    SccHHWx = cintersect(T, Yscc, HalfMaxScc); SccHHW = abs(diff(SccHHWx));
    
    %Extracting information on secundary peaks in the SCC ...
    [~, ~, SccXsecPeaks, SccYsecPeaks] = getPeaks(T, Yscc, 0, DomPer);

    %Reorganizing calculated data ...
    SXCC.lag        = T;
    SXCC.normco     = Yscc;
    [SXCC.max, idx] = max(Yscc);
    SXCC.rate       = max(YsccRate);
    SXCC.lagatmax   = T(idx);
    SXCC.maxsecpks  = SccYsecPeaks;
    SXCC.lagsecpks  = SccXsecPeaks;
    SXCC.hhw        = SccHHW;
    SXCC.hhwx       = SccHHWx;
    SXCC.halfmax    = HalfMaxScc;
    SXCC.fft.freq   = FFTscc.Freq;
    SXCC.fft.p      = FFTscc.Magn.P;   
    SXCC.fft.db     = FFTscc.Magn.dB;
    SXCC.fft.df     = FFTscc.DF;
    SXCC.fft.bw     = FFTscc.BW;
    
    DCC.lag          = [];
    DCC.normco       = [];
    DCC.max          = NaN;
    DCC.rate         = NaN;
    DCC.lagatmax     = NaN;
    DCC.maxsecpks    = [NaN, NaN];
    DCC.lagsecpks    = [NaN, NaN];
    DCC.fft.freq     = [];
    DCC.fft.p        = [];   
    DCC.fft.db       = [];
    DCC.fft.df       = NaN;
    DCC.fft.bw       = NaN;
    DCC.env.max      = NaN;
    DCC.env.lagatmax = NaN;
    DCC.env.hhw      = NaN;
    DCC.env.hhwx     = [NaN, NaN];
    DCC.env.halfmax  = NaN;
end

%----------------------------------------------------------------------------
function Y = ApplyNorm(Y, N, NormStr)

if (nargin == 2), NormStr = 'dries'; end

switch lower(NormStr)
case 'dries'
    if ~all(Y == 0)
        Y = Y/N.DriesNorm;
    else
        Y = ones(size(Y));
    end
case 'rate'
    Y = 1e3*Y/N.NF;
end

%----------------------------------------------------------------------------
function DF = DetermineCalcDF(ParamCalcDF, ThrCF, DifDF, SacDF)

if isnumeric(ParamCalcDF)
    if ~isnan(ParamCalcDF)
        DF = ParamCalcDF;
    elseif ~isnan(DifDF)
        DF = DifDF;
    elseif ~isnan(SacDF)
        DF = SacDF;
    else
        DF = ThrCF;
    end
elseif strcmpi(ParamCalcDF, 'cf')
    DF = ThrCF;
elseif strcmpi(ParamCalcDF, 'df') 
    if ~isnan(DifDF)
        DF = DifDF;
    else
        DF = SacDF;
    end
else
    DF = NaN;
end

%----------------------------------------------------------------------------
function PlotData(SXAC, DAC, MAC, SXCC, DCC, MCC, GBOR, Thr, RC, Info, StimParam, Param)

%Creating figure ...
FigHdl = figure('Name', sprintf('%s: %s', upper(mfilename), Info.capstr), ...
    'NumberTitle', 'off', ...
    'Units', 'normalized', ...
    'OuterPosition', [0 0.025 1 0.975], ... %Maximize figure (not in the MS Windows style!) ...
    'PaperType', 'A4', ...
    'PaperPositionMode', 'manual', ...
    'PaperUnits', 'normalized', ...
    'PaperPosition', [0.05 0.05 0.90 0.90], ...
    'PaperOrientation', 'portrait');

%Plot header ...
Str = { Info.hdrstr, sprintf('\\rm\\fontsize{9}Created by %s @ %s', upper(mfilename), datestr(now))};
AxHDR = axes('Position', [0.05 0.90 0.90 0.10], 'Visible', 'off');
text(0.5, 0.5, Str, 'Units', 'normalized', 'VerticalAlignment', 'middle', 'HorizontalAlignment', 'center', ...
    'FontWeight', 'bold', 'FontSize', 12);

%Plotting correlation functions ...
PlotAutoCorFnc([0.00, 0.0, 0.4, 0.8], SXAC(1), DAC(1), MAC(1), Param);
PlotAutoCorFnc([0.33, 0.0, 0.4, 0.8], SXAC(2), DAC(2), MAC(1), Param);
PlotCrossCorFnc([0.66, 0.0, 0.4, 0.8], SXCC, DCC, MCC, GBOR, Param);

%Display general information ...
PlotGenInfo([0.0, 0.73, 0.4, 0.3], Info, Thr, RC, GBOR, StimParam, Param);

%----------------------------------------------------------------------------
function PlotAutoCorFnc(ViewPort, SXAC, DAC, MAC, Param)

%General plot parameters ...
TitleFontSz  = 9;    %in points ...
LabelFontSz  = 7;    %in points ...
TckMrkFontSz = 7;    %in points ...
VerMargin    = 0.10; %in percent ...
HorMargin    = 0.15; %in percent ...
VerAxSpacing = 0.10; %in percent ...

Width     = ViewPort(3)*(1-2*HorMargin);
Height    = ViewPort(4)*(1-2*VerMargin-2*VerAxSpacing);
AxSpacing = ViewPort(4)*VerAxSpacing;
Origin    = ViewPort(1:2)+[Width*HorMargin, Height*VerMargin];

%Normalisation Coincidence count ...
NormStr = sprintf('Norm. Count\n(N_{rep}*(N_{rep}-1)*r^2*\\Delta\\tau*D)');

%Plotting SAC and if possible XAC ...
X = SXAC.lag; Y = SXAC.normco;
[MinX, MaxX, XTicks] = GetAxisLim('X', X, Param.corxrange, Param.corxstep);
Pos = [Origin(1), Origin(2)+3*AxSpacing+0.75*Height, Width, Height*0.25];
AxCOR = axes('Position', Pos, 'Box', 'off', 'TickDir', 'out', 'FontSize', TckMrkFontSz);
if (size(Y, 1) == 1)
    LnHdl = line(X, Y, 'LineStyle', '-', 'Marker', 'none', 'LineWidth', 2, 'Color', 'b');
    title('SAC', 'FontSize', TitleFontSz);
else
    LnHdl = line(X, Y, 'LineStyle', '-', 'Marker', 'none'); 
    set(LnHdl(1), 'LineWidth', 2, 'Color', 'b'); set(LnHdl(2), 'LineWidth', 0.5, 'Color', 'g');
    title('SAC and XAC', 'FontSize', TitleFontSz);
    LgHdl = legend({'SAC', 'XAC'});
    set(findobj(LgHdl, 'type', 'text'), 'FontSize', LabelFontSz);
end
TxtStr =  {sprintf('Max(SAC): %.2f', SXAC.max), sprintf('HHW(SAC): %.2fms', SXAC.hhw)};
xlabel('Delay(ms)', 'FontSize', LabelFontSz); ylabel(NormStr, 'Units', 'normalized', 'Position', [-0.11, 0.5, 0], 'FontSize', LabelFontSz);
set(AxCOR, 'XLim', [MinX MaxX], 'XTick', XTicks); 
YRange = get(AxCOR, 'YLim'); MinY = YRange(1); MaxY = YRange(2);
LnHdl = plotcintersect(SXAC.hhwx, SXAC.halfmax([1 1]), MinY);
set(LnHdl(1), 'LineStyle', '-', 'Color', [0 0 0], 'LineWidth', 0.5);
PlotVerZero(MinY, MaxY);
text(MinX, MaxY, TxtStr, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', 'FontSize', LabelFontSz);

%Plotting DIFCOR ...
Pos = [Origin(1), Origin(2)+2*AxSpacing+0.5*Height, Width, Height*0.25];
if isempty(DAC.lag)
    AxDIF = CreateEmptyAxis(Pos, LabelFontSz);
else
    X = DAC.lag; Y = DAC.normco; 
    AxDIF = axes('Position', Pos, 'Box', 'off', 'TickDir', 'out', 'FontSize', TckMrkFontSz);
    LnHdl = line(X, Y); 
    set(LnHdl(1), 'LineStyle', '-', 'Color', [1 0 0], 'LineWidth', 2);
    set(LnHdl(2), 'LineStyle', '-', 'Color', [0 0 0], 'LineWidth', 0.5);
    set(LnHdl(3), 'LineStyle', '-', 'Color', [0 0 0], 'LineWidth', 0.5);
    title('DIFCOR', 'FontSize', TitleFontSz);
    xlabel('Delay(ms)', 'FontSize', LabelFontSz); ylabel(NormStr, 'Units', 'normalized', 'Position', [-0.11, 0.5, 0], 'FontSize', LabelFontSz);
    set(AxDIF, 'XLim', [MinX MaxX], 'XTick', XTicks);  
    YRange = get(AxDIF, 'YLim'); MinY = YRange(1); MaxY = YRange(2);
    PlotVerZero(MinY, MaxY); PlotHorZero(MinX, MaxX);
    LnHdl = plotcintersect(DAC.env.hhwx, DAC.env.halfmax([1 1]), MinY);
    set(LnHdl(1), 'LineStyle', '-', 'Color', [0 0 0], 'LineWidth', 0.5);
    text(MinX, MaxY, {sprintf('Max(DIF): %.2f', DAC.max), sprintf('HHW(ENV): %.2fms', DAC.env.hhw)},'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', 'FontSize', LabelFontSz);
end

%Plotting SUMCOR ...
Pos = [Origin(1), Origin(2)+AxSpacing+0.25*Height, Width, Height*0.25];
if isempty(MAC.lag)
    AxDIF = CreateEmptyAxis(Pos, LabelFontSz);
else
    X = MAC.lag; Y = MAC.normco; 
    AxDIF = axes('Position', Pos, 'Box', 'off', 'TickDir', 'out', 'FontSize', TckMrkFontSz);
    LnHdl = line(X, Y); 
    set(LnHdl(1), 'LineStyle', '-', 'Color', [1 0 0], 'LineWidth', 2);
    set(LnHdl(2), 'LineStyle', '-', 'Color', [0 0 0], 'LineWidth', 0.5);
    set(LnHdl(3), 'LineStyle', '-', 'Color', [0 0 0], 'LineWidth', 0.5);
    title('SUMCOR', 'FontSize', TitleFontSz);
    xlabel('Delay(ms)', 'FontSize', LabelFontSz); ylabel(NormStr, 'Units', 'normalized', 'Position', [-0.11, 0.5, 0], 'FontSize', LabelFontSz);
    set(AxDIF, 'XLim', [MinX MaxX], 'XTick', XTicks);  
    YRange = get(AxDIF, 'YLim'); MinY = YRange(1); MaxY = YRange(2);
    PlotVerZero(MinY, MaxY); PlotHorZero(MinX, MaxX);
    LnHdl = plotcintersect(MAC.env.hhwx, MAC.env.halfmax([1 1]), MinY);
    set(LnHdl(1), 'LineStyle', '-', 'Color', [0 0 0], 'LineWidth', 0.5);
    text(MinX, MaxY, {sprintf('Max(SUM): %.2f', MAC.max), sprintf('HHW(ENV): %.2fms', MAC.env.hhw)},'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', 'FontSize', LabelFontSz);
end

%Plotting discrete fourier transform of DIFCOR and SUMCOR or SAC...
Pos = [Origin(1), Origin(2), Width, 0.25*Height];
if isempty(DAC.lag)
    FFT = SXAC.fft; TitleStr = 'DFT on SAC';
    X = FFT.freq;
    if strcmpi(Param.fftyunit, 'dB')
        Y = FFT.db;
        YLblStr = 'Power (dB, 10 log)'; 
    else
        Y = FFT.p;
        YLblStr = 'Power';
    end
    if ~isnan(FFT.df)
        Ord = floor(log10(FFT.df*2))-1;
        MinX = 0;
        MaxX = round(FFT.df*2*10^-Ord)*10^Ord;
        XTicks = 'auto';
    else
        [MinX, MaxX, XTicks] = GetAxisLim('X', X, Param.fftxrange, Param.fftxstep);
    end
    [MinY, MaxY] = GetAxisLim('Y', Y, Param.fftyrange);
    AxDFT = axes('Position', Pos, 'Box', 'off', 'TickDir', 'out', 'FontSize', TckMrkFontSz);
    LnHdl = line(X, Y);
    set(LnHdl(2), 'LineStyle', '-', 'Color', 'b', 'LineWidth', 0.5);
    set(LnHdl(1), 'LineStyle', ':', 'Color', 'k', 'LineWidth', 0.5);
    title(TitleStr, 'FontSize', TitleFontSz);
    xlabel('Freq(Hz)', 'FontSize', LabelFontSz); ylabel(YLblStr, 'FontSize', LabelFontSz);
    if ~ischar(XTicks)
        set(AxDFT, 'XLim', [MinX MaxX], 'YLim', [MinY MaxY], 'XTick', XTicks); 
    else
        set(AxDFT, 'XLim', [MinX MaxX], 'YLim', [MinY MaxY]);
    end
    YRange = get(AxDFT, 'YLim'); MinY = YRange(1); MaxY = YRange(2);
    line(FFT.df([1,1]), [MinY, MaxY], 'Color', [0 0 0], 'LineStyle', ':'); %Vertical line at dominant frequency ...
    text(MinX, MaxY, {sprintf('DomFreq: %.2fHz', FFT.df), sprintf('BandWidth: %.2fHz', FFT.bw)},'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', 'FontSize', LabelFontSz);
    legend({'Orig', 'RunAv'});
else
    FFT1 = DAC.fft; FFT2 = MAC.fft;
    TitleStr = 'DFT on DIFCOR and SUMCOR';
    
    X1 = FFT1.freq; X2 = FFT2.freq;
    if strcmpi(Param.fftyunit, 'dB')
        Y1 = FFT1.db; Y2 = FFT2.db;
        YLblStr = 'Power (dB, 10 log)'; 
    else
        Y1 = FFT1.p; Y2 = FFT2.p;
        YLblStr = 'Power';
    end
    if ~isnan(FFT1.df)
        Ord = floor(log10(FFT1.df*2))-1;
        MinX = 0;
        MaxX = round(FFT1.df*2*10^-Ord)*10^Ord;
        XTicks = 'auto';
    else
        [MinX1, MaxX1, ~] = GetAxisLim('X', X1, Param.fftxrange, Param.fftxstep);
        [MinX2, MaxX2, XTicks] = GetAxisLim('X', X2, Param.fftxrange, Param.fftxstep);
        MinX = min([MinX1, MinX2]); MaxX = max([MaxX1, MaxX2]);
    end
    [MinY1, MaxY1] = GetAxisLim('Y', Y1, Param.fftyrange);
    [MinY2, MaxY2] = GetAxisLim('Y', Y2, Param.fftyrange);
    MinY = min([MinY1, MinY2]); MaxY = max([MaxY1, MaxY2]);
    
    AxDFT = axes('Position', Pos, 'Box', 'off', 'TickDir', 'out', 'FontSize', TckMrkFontSz);
    LnHdl = line(X1, Y1);
    set(LnHdl(2), 'LineStyle', '-', 'Color', 'b', 'LineWidth', 0.5);
    set(LnHdl(1), 'LineStyle', ':', 'Color', 'k', 'LineWidth', 0.5);
    LnHdl = line(X2, Y2);
    set(LnHdl(2), 'LineStyle', '-', 'Color', 'r', 'LineWidth', 0.5);
    set(LnHdl(1), 'LineStyle', ':', 'Color', 'k', 'LineWidth', 0.5);
    
    title(TitleStr, 'FontSize', TitleFontSz);
    xlabel('Freq(Hz)', 'FontSize', LabelFontSz); ylabel(YLblStr, 'FontSize', LabelFontSz);
    if ~ischar(XTicks)
        set(AxDFT, 'XLim', [MinX MaxX], 'YLim', [MinY MaxY], 'XTick', XTicks); 
    else
        set(AxDFT, 'XLim', [MinX MaxX], 'YLim', [MinY MaxY]);
    end
    YRange = get(AxDFT, 'YLim'); MinY = YRange(1); MaxY = YRange(2);
    line(FFT1.df([1,1]), [MinY, MaxY], 'Color', [0 0 0], 'LineStyle', ':'); %Vertical line at dominant frequency ...
    line(FFT2.df([1,1]), [MinY, MaxY], 'Color', [0 0 0], 'LineStyle', ':'); %Vertical line at dominant frequency ...
    text(MinX, MaxY, {sprintf('DomFreq Dif: %.2fHz, DomFreq Sum: %.2fHz', FFT1.df, FFT2.df), ... 
        sprintf('BandWidth Dif: %.2fHz, BandWidth Sum: %.2fHz', FFT1.bw, FFT2.bw)}, ...
        'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', 'FontSize', LabelFontSz);
    legend({'Orig Dif', 'RunAv Dif', 'Orig Sum', 'RunAv Sum'});
end

%----------------------------------------------------------------------------
function PlotCrossCorFnc(ViewPort, SXCC, DCC, MCC, GBOR, Param)

%General plot parameters ...
TitleFontSz  = 9;    %in points ...
LabelFontSz  = 7;    %in points ...
TckMrkFontSz = 7;    %in points ...
VerMargin    = 0.10; %in percent ...
HorMargin    = 0.15; %in percent ...
VerAxSpacing = 0.10; %in percent ...

Width     = ViewPort(3)*(1-2*HorMargin);
Height    = ViewPort(4)*(1-2*VerMargin-2*VerAxSpacing);
AxSpacing = ViewPort(4)*VerAxSpacing;
Origin    = ViewPort(1:2)+[Width*HorMargin, Height*VerMargin];

%Normalisation Coincidence count ...
NormStr = sprintf('Norm. Count\n(N_{rep}*(N_{rep}-1)*r^2*\\Delta\\tau*D)');

%Plotting SCC and if possible XCC ...
X = SXCC.lag; Y = SXCC.normco;
[MinX, MaxX, XTicks] = GetAxisLim('X', X, Param.corxrange, Param.corxstep);
Pos = [Origin(1), Origin(2)+3*AxSpacing+0.75*Height, Width, Height*0.25];
AxCOR = axes('Position', Pos, 'Box', 'off', 'TickDir', 'out', 'FontSize', TckMrkFontSz);
if (size(Y, 1) == 1)
    LnHdl = line(X, Y, 'LineStyle', '-', 'Marker', 'none', 'LineWidth', 2, 'Color', 'b');
    title('SCC', 'FontSize', TitleFontSz);
else
    LnHdl = line(X, Y, 'LineStyle', '-', 'Marker', 'none'); 
    set(LnHdl(1), 'LineWidth', 2, 'Color', 'b'); set(LnHdl(2), 'LineWidth', 0.5, 'Color', 'g');
    title('SCC and XCC', 'FontSize', TitleFontSz);
    LgHdl = legend({'SCC', 'XCC'});
    set(findobj(LgHdl, 'type', 'text'), 'FontSize', LabelFontSz);
end
TxtStr =  {sprintf('Max(SCC): %.2f @ %.2fms', SXCC.max, SXCC.lagatmax), sprintf('HHW(SCC): %.2fms', SXCC.hhw)};
xlabel('Delay(ms)', 'FontSize', LabelFontSz); ylabel(NormStr, 'Units', 'normalized', 'Position', [-0.11, 0.5, 0], 'FontSize', LabelFontSz);
set(AxCOR, 'XLim', [MinX MaxX], 'XTick', XTicks); 
YRange = get(AxCOR, 'YLim'); MinY = YRange(1); MaxY = YRange(2);
LnHdl = plotcintersect(SXCC.hhwx, SXCC.halfmax([1 1]), MinY);
set(LnHdl(1), 'LineStyle', '-', 'Color', [0 0 0], 'LineWidth', 0.5);
line(SXCC.lagatmax([1 1]), [MinY, MaxY], 'Color', [0 0 0], 'LineStyle', ':');
line(SXCC.lagatmax, SXCC.max, 'Color', [0 0 0], 'LineStyle', 'none', 'marker', '.');
line(SXCC.lagsecpks([1 1]), [MinY, MaxY], 'Color', [0 0 0], 'LineStyle', ':');
line(SXCC.lagsecpks([2 2]), [MinY, MaxY], 'Color', [0 0 0], 'LineStyle', ':');
line(SXCC.lagsecpks, SXCC.maxsecpks, 'Color', [0 0 0], 'LineStyle', 'none', 'marker', '.');
PlotVerZero(MinY, MaxY);
text(MinX, MaxY, TxtStr, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', 'FontSize', LabelFontSz);

%Plotting DIFCOR ...
Pos = [Origin(1), Origin(2)+2*AxSpacing+0.5*Height, Width, Height*0.25];
if isempty(DCC.lag)
    AxDIF = CreateEmptyAxis(Pos, LabelFontSz);
else
    X = DCC.lag; Y = DCC.normco; 
    AxDIF = axes('Position', Pos, 'Box', 'off', 'TickDir', 'out', 'FontSize', TckMrkFontSz);
    LnHdl = line(X, Y); 
	if Param.gaborfit && ~isempty(GBOR) && isempty(GBOR.err)
		LnHdl(end+1) = line(GBOR.x, GBOR.y, 'LineStyle', '--');
		LnHdl(end+1) = line(GBOR.x, GBOR.env, 'LineStyle', '--');
		LnHdl(end+1) = line(GBOR.x, -GBOR.env, 'LineStyle', '--');
	end
    set(LnHdl(1), 'LineStyle', '-', 'Color', [1 0 0], 'LineWidth', 2);
    set(LnHdl(2), 'LineStyle', '-', 'Color', [0 0 0], 'LineWidth', 0.5);
    set(LnHdl(3), 'LineStyle', '-', 'Color', [0 0 0], 'LineWidth', 0.5);
    title('DIFCOR', 'FontSize', TitleFontSz);
    xlabel('Delay(ms)', 'FontSize', LabelFontSz); 
	ylabel(NormStr, 'Units', 'normalized', 'Position', [-0.11, 0.5, 0], 'FontSize', LabelFontSz);
    set(AxDIF, 'XLim', [MinX MaxX], 'XTick', XTicks);  
    YRange = get(AxDIF, 'YLim');
    MinY = YRange(1); MaxY = YRange(2);
    PlotVerZero(MinY, MaxY);
    PlotHorZero(MinX, MaxX);
    LnHdl = plotcintersect(DCC.env.hhwx, DCC.env.halfmax([1 1]), MinY);
    set(LnHdl(1), 'LineStyle', '-', 'Color', [0 0 0], 'LineWidth', 0.5);
    line(DCC.lagatmax, DCC.max, 'Color', [0 0 0], 'LineStyle', 'none', 'marker', '.');
    line(DCC.lagatmax([1 1]), [MinY, MaxY], 'Color', [0 0 0], 'LineStyle', ':');
    line(DCC.lagsecpks([1 1]), [MinY, MaxY], 'Color', [0 0 0], 'LineStyle', ':');
    line(DCC.lagsecpks([2 2]), [MinY, MaxY], 'Color', [0 0 0], 'LineStyle', ':');
    line(DCC.lagsecpks, DCC.maxsecpks, 'Color', [0 0 0], 'LineStyle', 'none', 'marker', '.');
    text(MinX, MaxY, {sprintf('Max(DIF): %.2f @ %.2fms', DCC.max, DCC.lagatmax), ...
            sprintf('SecPks(DIF): @ %.2fms & %.2fms', DCC.lagsecpks), ...
            sprintf('Max(ENV):  %.2f @ %.2fms', DCC.env.max, DCC.env.lagatmax), ...
            sprintf('HHW(ENV): %.2fms', DCC.env.hhw)},'HorizontalAlignment', ...
            'left', 'VerticalAlignment', 'top', 'FontSize', LabelFontSz);
end

%Plotting SUMCOR ...
Pos = [Origin(1), Origin(2)+AxSpacing+0.25*Height, Width, Height*0.25];
if isempty(MCC.lag)
    AxDIF = CreateEmptyAxis(Pos, LabelFontSz);
else
    X = MCC.lag; Y = MCC.normco; 
    AxDIF = axes('Position', Pos, 'Box', 'off', 'TickDir', 'out', 'FontSize', TckMrkFontSz);
    LnHdl = line(X, Y); 
	if Param.gaborfit && ~isempty(GBOR) && isempty(GBOR.err)
		LnHdl(end+1) = line(GBOR.x, GBOR.y, 'LineStyle', '--');
		LnHdl(end+1) = line(GBOR.x, GBOR.env, 'LineStyle', '--');
		LnHdl(end+1) = line(GBOR.x, -GBOR.env, 'LineStyle', '--');
	end
    set(LnHdl(1), 'LineStyle', '-', 'Color', [1 0 0], 'LineWidth', 2);
    set(LnHdl(2), 'LineStyle', '-', 'Color', [0 0 0], 'LineWidth', 0.5);
    set(LnHdl(3), 'LineStyle', '-', 'Color', [0 0 0], 'LineWidth', 0.5);
    title('SUMCOR', 'FontSize', TitleFontSz);
    xlabel('Delay(ms)', 'FontSize', LabelFontSz); 
	ylabel(NormStr, 'Units', 'normalized', 'Position', [-0.11, 0.5, 0], 'FontSize', LabelFontSz);
    set(AxDIF, 'XLim', [MinX MaxX], 'XTick', XTicks);  
    YRange = get(AxDIF, 'YLim');
    MinY = YRange(1); MaxY = YRange(2);
    PlotVerZero(MinY, MaxY);
    PlotHorZero(MinX, MaxX);
    LnHdl = plotcintersect(MCC.env.hhwx, MCC.env.halfmax([1 1]), MinY);
    set(LnHdl(1), 'LineStyle', '-', 'Color', [0 0 0], 'LineWidth', 0.5);
    line(MCC.lagatmax, MCC.max, 'Color', [0 0 0], 'LineStyle', 'none', 'marker', '.');
    line(MCC.lagatmax([1 1]), [MinY, MaxY], 'Color', [0 0 0], 'LineStyle', ':');
    line(MCC.lagsecpks([1 1]), [MinY, MaxY], 'Color', [0 0 0], 'LineStyle', ':');
    line(MCC.lagsecpks([2 2]), [MinY, MaxY], 'Color', [0 0 0], 'LineStyle', ':');
    line(MCC.lagsecpks, MCC.maxsecpks, 'Color', [0 0 0], 'LineStyle', 'none', 'marker', '.');
    text(MinX, MaxY, {sprintf('Max(SUM): %.2f @ %.2fms', DCC.max, DCC.lagatmax), ...
            sprintf('SecPks(SUM): @ %.2fms & %.2fms', DCC.lagsecpks), ...
            sprintf('Max(ENV):  %.2f @ %.2fms', DCC.env.max, DCC.env.lagatmax), ...
            sprintf('HHW(ENV): %.2fms', DCC.env.hhw)},'HorizontalAlignment', ...
            'left', 'VerticalAlignment', 'top', 'FontSize', LabelFontSz);
end

%Plotting discrete fourier transform of DIFCOR and SUMCOR or SAC...
Pos = [Origin(1), Origin(2), Width, 0.25*Height];
if isempty(DCC.lag)
    FFT = SXCC.fft; TitleStr = 'DFT on SCC';
    X = FFT.freq;
    if strcmpi(Param.fftyunit, 'dB')
        Y = FFT.db;
        YLblStr = 'Power (dB, 10 log)'; 
    else
        Y = FFT.p;
        YLblStr = 'Power';
    end
    if ~isnan(FFT.df)
        Ord = floor(log10(FFT.df*2))-1;
        MinX = 0;
        MaxX = round(FFT.df*2*10^-Ord)*10^Ord;
        XTicks = 'auto';
    else
        [MinX, MaxX, XTicks] = GetAxisLim('X', X, Param.fftxrange, Param.fftxstep);
    end
    [MinY, MaxY] = GetAxisLim('Y', Y, Param.fftyrange);
    AxDFT = axes('Position', Pos, 'Box', 'off', 'TickDir', 'out', 'FontSize', TckMrkFontSz);
    LnHdl = line(X, Y);
    set(LnHdl(2), 'LineStyle', '-', 'Color', 'b', 'LineWidth', 0.5);
    set(LnHdl(1), 'LineStyle', ':', 'Color', 'k', 'LineWidth', 0.5);
    title(TitleStr, 'FontSize', TitleFontSz);
    xlabel('Freq(Hz)', 'FontSize', LabelFontSz); ylabel(YLblStr, 'FontSize', LabelFontSz);
    if ~ischar(XTicks)
        set(AxDFT, 'XLim', [MinX MaxX], 'YLim', [MinY MaxY], 'XTick', XTicks); 
    else
        set(AxDFT, 'XLim', [MinX MaxX], 'YLim', [MinY MaxY]);
    end
    YRange = get(AxDFT, 'YLim'); MinY = YRange(1); MaxY = YRange(2);
    line(FFT.df([1,1]), [MinY, MaxY], 'Color', [0 0 0], 'LineStyle', ':'); %Vertical line at dominant frequency ...
    text(MinX, MaxY, {sprintf('DomFreq: %.2fHz', FFT.df), sprintf('BandWidth: %.2fHz', FFT.bw)},'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', 'FontSize', LabelFontSz);
    legend({'Orig', 'RunAv'});
else
    FFT1 = DCC.fft; FFT2 = MCC.fft;
    TitleStr = 'DFT on DIFCOR and SUMCOR';
    
    X1 = FFT1.freq; X2 = FFT2.freq;
    if strcmpi(Param.fftyunit, 'dB')
        Y1 = FFT1.db; Y2 = FFT2.db;
        YLblStr = 'Power (dB, 10 log)'; 
    else
        Y1 = FFT1.p; Y2 = FFT2.p;
        YLblStr = 'Power';
    end
    if ~isnan(FFT1.df)
        Ord = floor(log10(FFT1.df*2))-1;
        MinX = 0;
        MaxX = round(FFT1.df*2*10^-Ord)*10^Ord;
        XTicks = 'auto';
    else
        [MinX1, MaxX1, ~] = GetAxisLim('X', X1, Param.fftxrange, Param.fftxstep);
        [MinX2, MaxX2, XTicks] = GetAxisLim('X', X2, Param.fftxrange, Param.fftxstep);
        MinX = min([MinX1, MinX2]); MaxX = max([MaxX1, MaxX2]);
    end
    [MinY1, MaxY1] = GetAxisLim('Y', Y1, Param.fftyrange);
    [MinY2, MaxY2] = GetAxisLim('Y', Y2, Param.fftyrange);
    MinY = min([MinY1, MinY2]); MaxY = max([MaxY1, MaxY2]);
    
    AxDFT = axes('Position', Pos, 'Box', 'off', 'TickDir', 'out', 'FontSize', TckMrkFontSz);
    LnHdl = line(X1, Y1);
    set(LnHdl(2), 'LineStyle', '-', 'Color', 'b', 'LineWidth', 0.5);
    set(LnHdl(1), 'LineStyle', ':', 'Color', 'k', 'LineWidth', 0.5);
    LnHdl = line(X2, Y2);
    set(LnHdl(2), 'LineStyle', '-', 'Color', 'r', 'LineWidth', 0.5);
    set(LnHdl(1), 'LineStyle', ':', 'Color', 'k', 'LineWidth', 0.5);
    
    title(TitleStr, 'FontSize', TitleFontSz);
    xlabel('Freq(Hz)', 'FontSize', LabelFontSz); ylabel(YLblStr, 'FontSize', LabelFontSz);
    if ~ischar(XTicks)
        set(AxDFT, 'XLim', [MinX MaxX], 'YLim', [MinY MaxY], 'XTick', XTicks); 
    else
        set(AxDFT, 'XLim', [MinX MaxX], 'YLim', [MinY MaxY]);
    end
    YRange = get(AxDFT, 'YLim'); MinY = YRange(1); MaxY = YRange(2);
    line(FFT1.df([1,1]), [MinY, MaxY], 'Color', [0 0 0], 'LineStyle', ':'); %Vertical line at dominant frequency ...
    line(FFT2.df([1,1]), [MinY, MaxY], 'Color', [0 0 0], 'LineStyle', ':'); %Vertical line at dominant frequency ...
    text(MinX, MaxY, {sprintf('DomFreq Dif: %.2fHz, DomFreq Sum: %.2fHz', FFT1.df, FFT2.df), ... 
        sprintf('BandWidth Dif: %.2fHz, BandWidth Sum: %.2fHz', FFT1.bw, FFT2.bw)}, ...
        'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', 'FontSize', LabelFontSz);
    legend({'Orig Dif', 'RunAv Dif', 'Orig Sum', 'RunAv Sum'});
end


%----------------------------------------------------------------------------
function PlotGenInfo(ViewPort, Info, Thr, RC, GBOR, StimParam, Param)

AxINF = axes('Position', ViewPort, 'Visible', 'off');
text(0.35, 0.85, Info.hdrstr, 'Units', 'normalized', 'VerticalAlignment', 'middle', 'HorizontalAlignment', 'center', ...
    'FontWeight', 'demi', 'FontSize', 9);

text(0.05, 0.75, char(Thr(1).str, '', RC(1).str), 'Units', 'normalized', 'VerticalAlignment', 'middle', 'HorizontalAlignment', 'left', ...
    'FontWeight', 'normal', 'FontSize', 8);
text(0.35, 0.75, char(Thr(2).str, '', RC(2).str), 'Units', 'normalized', 'VerticalAlignment', 'middle', 'HorizontalAlignment', 'left', ...
    'FontWeight', 'normal', 'FontSize', 8);

text(0.05, 0.60, {'Stimulus parameters:', 'A^+_1, A^-_1, A^+_2, A^-_2'},'Units', 'normalized', 'VerticalAlignment', 'middle', 'HorizontalAlignment', 'left', ...
    'FontWeight', 'demi', 'FontSize', 9);
text(0.05, 0.45, StimParam.str, 'Units', 'normalized', 'VerticalAlignment', 'middle', 'HorizontalAlignment', 'left', ...
    'FontWeight', 'normal', 'FontSize', 8);

text(0.35, 0.60, 'Calculation parameters:','Units', 'normalized', 'VerticalAlignment', 'middle', 'HorizontalAlignment', 'left', ...
    'FontWeight', 'demi', 'FontSize', 9);
text(0.35, 0.40, Param.str, 'Units', 'normalized', 'VerticalAlignment', 'middle', 'HorizontalAlignment', 'left', ...
    'FontWeight', 'normal', 'FontSize', 8);

if Param.gaborfit
	gStr = AssembleGABORStr(GBOR);
	text(0.65, 0.50, gStr, 'Units', 'normalized', 'VerticalAlignment', 'middle', 'HorizontalAlignment', 'left', ...
    'FontWeight', 'normal', 'FontSize', 8);
end


%----------------------------------------------------------------------------
function [MinVal, MaxVal, Ticks] = GetAxisLim(AxisType, Values, Range, Step)

if strcmpi(AxisType, 'x') %Abcissa ...
    Margin = 0.00;
    
    if isinf(Range(1)), MinVal = min(Values(:))*(1-Margin); else, MinVal = Range(1); end
    if isinf(Range(2)), MaxVal = max(Values(:))*(1+Margin); else, MaxVal = Range(2); end
    Ticks = MinVal:Step:MaxVal;
else %Ordinate ...    
    Margin = 0.05;
    
    if isinf(Range(1)), MinVal = min([0; Values(:)])*(1-Margin); else, MinVal = Range(1); end
    if isinf(Range(2)), MaxVal = max(Values(:))*(1+Margin); else, MaxVal = Range(2); end
end

%----------------------------------------------------------------------------
function PlotVerZero(MinY, MaxY)

line([0, 0], [MinY, MaxY], 'Color', [0 0 0], 'LineStyle', ':');

%----------------------------------------------------------------------------
function PlotHorZero(MinX, MaxX)

line([MinX, MaxX], [0, 0], 'Color', [0 0 0], 'LineStyle', ':');

%----------------------------------------------------------------------------
function AxHdl = CreateEmptyAxis(Pos, FontSz)

%Create axis object ... 
AxHdl = axes('Position', Pos, 'Box', 'on', 'Color', [0.8 0.8 0.8], 'Units', 'normalized', ...
    'YTick', [], 'YTickLabel', [], 'XTick', [], 'XTickLabel', []);

%Create Text object ...
TxtHdl = text(0.5, 0.5, sprintf('No plot to be\ngenerated'));
set(TxtHdl, 'VerticalAlignment', 'middle', ...
             'HorizontalAlignment', 'center', ...
             'Color', 'k', ...
             'FontAngle', 'normal', ...
             'FontSize', FontSz, ...
             'FontWeight', 'normal');
         
%----------------------------------------------------------------------------        
function GABOR = calcGABOR(DCC, Param, Freq)

Period = 1000/Freq;
if isempty(Param.fitrange)
	Param.fitrange = Param.corxrange * 1000;
end

%Initial values for gabor fit
gaborParm.ampl = DCC.env.max;
maxEst = DCC.env.lagatmax;
gaborParm.max = [ maxEst-(1000/DCC.fft.df), maxEst,  maxEst+(1000/DCC.fft.df) ];
gaborParm.width = DCC.env.hhw;
gaborParm.freq = [DCC.fft.df/2 DCC.fft.df DCC.fft.df*2];
gaborParm.freq = gaborParm.freq/1000;
gaborParm.ph = ((DCC.lagatmax)*(DCC.fft.df/1000))*2*pi *-1;

%Run gabor fit
[X, Y, Env, Constants] = gaborfit(DCC.lag, DCC.normco(1,:), Param.fitrange/1000, Param.samplerate*1000, gaborParm); 

%get peaks
[BestITD, Max, SecPeaks] = getpeaks(X, Y, 0, Period);
BestITDc = BestITD/Period;
Peaks = [SecPeaks, BestITD];
[~, idx] = min(abs(Peaks));
ZeroPeak = Peaks(idx);

[EnvPeak, EnvMax] = getmaxloc(X, Env); %EnvPeak = EnvPeak/1000;
EnvPeakc = EnvPeak/Period;
InterSect = cintersect(X, Env, EnvMax/2); 
HHW  = diff(InterSect);
HHWc = HHW/Period;

CarrFreq = Constants.Freq * 1e3;
AccFrac = Constants.AccFraction;

GABOR = RecLowerFields(CollectInStruct(X, Y, Env, Constants, Max, BestITD, BestITDc, ...
        SecPeaks, ZeroPeak, EnvMax, EnvPeak, EnvPeakc, InterSect, HHW, HHWc, CarrFreq, AccFrac));
	
%----------------------------------------------------------------------------
function S = RecLowerFields(S)
FNames  = fieldnames(S);
NFields = length(FNames);
for n = 1:NFields
	Val = getfield(S, FNames{n});
	S = rmfield(S, FNames{n});
	if isstruct(Val), S = setfield(S, lower(FNames{n}), RecLowerFields(Val));  %Recursive ...
	else, S = setfield(S, lower(FNames{n}), Val); end
end
		
%----------------------------------------------------------------------------
function Str = AssembleGABORStr(GABOR)
if ~isempty(GABOR.err)
	Str = sprintf('\\bfGABOR FAILED:\\rm\n %s', GABOR.err);
	return
end

Str = {sprintf('\\bfGABOR\\rm'); ...
        sprintf('\\itGeneral\\rm'); ...
        sprintf('Max = %s', Param2Str(GABOR.max, 'spk/sec', 0)); ...
        sprintf('ITD = %s/%s', Param2Str(GABOR.bestitd, 'ms', 2), Param2Str(GABOR.bestitdc, 'cyc', 2)); ...
        sprintf('HHW = %s/%s', Param2Str(GABOR.hhw, 'ms', 2), Param2Str(GABOR.hhwc, 'cyc', 2)); ...
        sprintf('Peaks = %.2f/%.2f/%.2fms', [GABOR.secpeaks(1), GABOR.zeropeak, GABOR.secpeaks(2)]); ...
        sprintf('\\itConstants\\rm'); ...
        sprintf('Amp = %s', Param2Str(GABOR.constants.a, 'spk/sec', 0)); ...
        sprintf('EnvMax = %s', Param2Str(GABOR.constants.envmax, 'ms', 2)); ...
        sprintf('EnvWidth = %s', Param2Str(GABOR.constants.envwidth, 'ms', 2)); ...
        sprintf('Freq = %s', Param2Str(GABOR.constants.freq*1e3, 'kHz', 2)); ...
        sprintf('Phase = %s', Param2Str(GABOR.constants.ph/2/pi, 'cyc', 2)); ...
        sprintf('AccFrac = %s', Param2Str(GABOR.constants.accfraction*100, '%', 0)); ... 
        sprintf('\\itEnveloppe\\rm'); ...
        sprintf('Max = %s', Param2Str(GABOR.envmax, 'spk/sec', 0)); ...
        sprintf('Peak = %s/%s', Param2Str(GABOR.envpeak, 'ms', 2), Param2Str(GABOR.envpeakc, 'cyc', 2))};