function T = EvalISI(ds,varargin)
% EVALISI - Evaluate inter spike intervals
%   EVALISI(ds) evaluates the inter spike intervals in dataset ds.
%
%   More options can be given and will be passed to CalcISIH, see help
%   CalcISIH and EvalISI('factory') for detailed information.
%
%   E.g.  ds=dataset('A0242','86-8-SPL');
%         EvalISI(ds)

%% ---------------- CHANGELOG -----------------------

% Sunday Feb 8 2015 Tina
%   bugfix in the plot (yes/no) option. "strcmpi(Param.plot, 'yes')"
%   Increased output of T tab, including minimum ISI and mean or average ISI
%   values.
%   bugfix, Output T now compatible to view with "structview" function.


%% Params
DefParam.isubseqs    = 'all';        %The subsequences included in the analysis.
                                     %By default all subsequences are included ('all') ...
DefParam.ireps       = 'all';        %The repetitions for the subsequences included 
                                     %in the analysis. By default all repetitions
                                     %are included ('all'). For multiple subsequences
                                     %this can be a cell-array with each element a
                                     %numerical vector representing the repetition
                                     %numbers to include for a particular subsequence ...
DefParam.anwin       = 'burstdur';   %Analysis window in ms. This must be
                                     %a vector with an even number of elements.
                                     %Each pair of this vector designates a time-
                                     %interval included in the calculation. By default
                                     %this is from 0 to stimulus duration ('burstdur') ...
DefParam.viewport    = [0 10];       %Viewport, i.e. the actual view on the histogram (in ms) ...
DefParam.nbin        = 64;           %Number of bins in the requested viewport ...
DefParam.runav       = 0;            %Number of bins used in smoothing the histogram ...
                                     %Smoothing does not change the extracted calculation
                                     %parameters ...
DefParam.minisi      = 0;            %Minimum interspike interval in ms ...
DefParam.timesubtr   = 0;            %Subtraction of a constant time in ms from
                                     %all spiketimes ...
DefParam.plot        = 'yes';        %Produce a plot or not.('yes'/'no')

if (nargin == 1) && ischar(ds) && strcmpi(ds, 'factory')
    disp(DefParam);
    return;
end
if ~isa(ds,'dataset')
    error('First argument must be a dataset');
end
if isTHRdata(ds), 
    error('ISI histogram cannot be plotted for current dataset stimulus type'); 
end
ds = FillDataset(ds);

%Checking additional property-list
Param = checkproplist(DefParam, varargin{:});

%% Calculation
calcData = CalcISIH(ds,'isubseqs',Param.isubseqs,'ireps',Param.ireps,...
    'anwin',Param.anwin,'viewport',Param.viewport,'nbin',Param.nbin,...
    'runav',Param.runav,'minisi',Param.minisi,'timesubtr',Param.timesubtr);

%% Plot
nsub = length(calcData.hist);
bincenters = cell(1, nsub);
rates = cell(1, nsub);
for n = 1:nsub
    bincenters{n} = calcData.hist(n).bincenters;
    rates{n} = calcData.hist(n).rate;
end

% GridPlot(bincenters, rates, ds, 'xlabel', {'ISI(ms)'}, ...
%     'ylabel', {'rate(#Intervals/sec)'}, 'mfileName', mfilename, ...
%     'plotTypeHdl', @histPlotObject);
%below modified by Tina 18/2/2014
if strcmpi(Param.plot, 'yes')
GridPlot(bincenters, rates, ds, 'xlabel', {'ISI(ms)'}, ...
    'ylabel', {'rate(#Intervals/sec)'}, 'mfileName', mfilename, ...
    'plotTypeHdl', @histPlotObject,'nsub',nsub);
end

%% output struct
T.ds.filename = ds.ID.Experiment.ID.Name;
T.ds.icell = ds.ID.iCell;
T.ds.iseq = ds.ID.iDataset;
T.ds.seqid = [num2str(ds.ID.iCell) '-' num2str(ds.ID.iRecOfCell) '-' ds.StimType];
T.ds.dur = ds.Stim.Presentation.TotDur;
T.fcar = ds.Stim.Fcar;
T.nsub=ds.Stim.Presentation.Ncond;
T.Hisbincenters=bincenters;
T.Hisrates=rates;

T.createdby = mfilename;
T.params=calcData.param;

for i=1:length(calcData.hist)
    %all of the minisi from each spl
    minisiSPL(i)=calcData.hist(i).minisi;
end
% %By Abel: Add stimulus info
% %data may be sorted, which is reflected in the order of isubsequence.
% idx = calcData.param.isubseqs;
% 
% stimParam.spl = sortByIdx_(GetSPL(ds), idx); %stim spls
% stimParam.spl = stimParam.spl(:).';
% stimParam.fcar = sortByIdx_(ds.fcar, idx); %stim freq
% T.stim = stimParam;
% 
% T.minisiSPL=sortByIdx_(minisiSPL, idx); %how to get the peak value?
T.Avminisi = mean(minisiSPL);


% 
% function sortedValues = sortByIdx_(values, idx)
% if length(values) == 1
% 	sortedValues = values;
% 	return
% end
% 
% sortedValues = values(idx);



