function [H,G] = FOISI(D,figh,P)
% dataset/FOISI - (first order) inter spike interval histogram of a dataset
%    FOISI(D) calculates inter spike interval histograms for the spike times in
%    dataset D.
%
%    FOISI(D,figh) uses figure handle figh for plotting
%    (default = [] -> gcf). 
%
%    FOISI(D,figh,P) uses parameters P for displaying the FOISI.
%    P is typically a dataviewparam object or a valid 2nd input argument to
%    the dataviewparam constructor method, such as a parameter filename.
%
%    FOISI is a standard "dataviewer", meaning that it may serve as
%    viewer for online data analysis during data collection. In addition,
%    the plot generated by all dataviewers allow an interactive change of
%    analysis parameter view the Params|Edit pullodwn menu (Ctr-Q).
%    For details on dataviewers, see dataviewparam.
%
%    See also dataviewparam, dataset/enableparamedit.

% Handle the special case of parameter queries. 
% Do this immediately to avoid endless recursion with dataviewparam.
if isvoid(D) && isequal('params', figh),
    [H,G] = local_ParamGUI;
    return;
end

% Should we open a new figure or use an existing one?
if nargin<2 || isempty(figh),
    open_new = isempty(get(0,'CurrentFigure'));
    figh = gcf; 
else
    open_new = isSingleHandle(figh);
end

% Parameters
if nargin<3, P = []; end
if isempty(P), % use default paremeter set for this dataviewer
    P = dataviewparam(mfilename); 
end

% delegate the real work to local fcn
H = local_FOISI(D, figh, open_new, P);

% enable parameter editing when viewing offline
if isSingleHandle(figh, 'figure'), enableparamedit(D, P, figh); end;



%============================================================
%============================================================
function data_struct = local_FOISI(D, figh, open_new, P);
% the real work for the inter spike interval histogram
if isSingleHandle(figh, 'figure')
    figure(figh); clf; ah = gca;
    if open_new, placefig(figh, mfilename, D.Stim.GUIname); end % restore previous size 
else
    ah = axes('parent', figh);
end

% Check varied stimulus Params
Pres = D.Stim.Presentation;
P = struct(P); P = P.Param;
isortPlot = P.iCond(P.iCond<=Pres.Ncond); % limit to actual Ncond
if isortPlot==0, isortPlot = 1:Pres.Ncond; end;
Ncond = numel(isortPlot);
AW = P.Anwin;
if isequal(P.scaling,'constrained'), doScale = 1; else doScale = 0; end;
% XXXXXX

% prepare plot
Clab = cellify(CondLabel(D));
[axh, Lh, Bh] = plotpanes(Ncond+1, 0, figh);

% get sorted spikes
Chan = 1; % digital input
TC = spiketimes(D, Chan, 'no-unwarp');

H = zeros(Ncond, P.Nbin);
isortPlot=isortPlot(:).';
for i=1:Ncond
    icond = isortPlot(i);
    BurstDur = max(burstdur(D,icond));
    if isequal('burstdur', AW), aw = [0 BurstDur]; else aw = AW; end; 
    if isequal('burstdur', P.viewport), P.viewport = [0 5]; end;
    BinEdges = linspace(P.viewport(1), P.viewport(2), P.Nbin+1); % histogram bins
    BinWidth = BinEdges(2)-BinEdges(1);
    BinCenters = BinEdges(1:(end-1))+(BinWidth/2);
    DiffSpt = [];
    for n = 1:Pres.Nrep % spike time differences must be calculated with a for loop
        spt = TC{icond,n}; % spike times of condition icond and repetition n
        spt = AnWin(spt, aw); % apply analysis window
        DiffSpt = [DiffSpt diff(spt)];
    end
    if isempty(DiffSpt)
        DiffSpt = Inf; 
    end % to avoid crash dump but still use histc
    N = histc(DiffSpt, BinEdges);
    N(end) = []; % remove last garbage bin
    H(i,:) = N;
    Mean = sum(BinCenters.*N)/sum(N); % mean interval in ms
    Rate = 1e3/Mean; % mean firing rate [/s]
%     Rstr = ['rate = ' num2str(Rate) ' spikes/s']; replaced (Jan 2018)
    Rstr = ['rate = ' num2str(round(Rate, 2)) ' spikes/s']; % round on two decimals (Jan 2018)
    h = axh(i); % current axes handle
    % axes(h); % slow!!!!
    bar(h, BinCenters, N, 'histc');
    xlim(h, [P.viewport(1) P.viewport(2)])
    title(h,{Clab{icond},Rstr},'FontSize', 8); % Rstr added to te title, more readable (Jan 2018)
%     title(h, Clab{icond});
    set(gcf,'CurrentAxes',h);
%     text(0.1, 0.1, Rstr, 'units', 'normalized', 'color', 'r', 'fontsize', 12 , 'interpreter', 'latex'); % This ifno is put in the title (Jan 2018)

    data_struct.BurstDur = BurstDur;
    data_struct.aw = aw;
    data_struct.BinEdges(icond,:) = BinEdges;
    data_struct.BinWidth(icond,:) = BinWidth;
    data_struct.BinCenters(icond,:) = BinCenters;
    data_struct.DiffSpt{icond} = DiffSpt;
    data_struct.HistogramCount{icond} = N;
    data_struct.MeanInterval(icond) = Mean;
    data_struct.Rate(icond) = Rate;
    data_struct.Rstr{icond} = Rstr;
    data_struct.xlim = [P.viewport(1) P.viewport(2)];
    data_struct.title{icond} = Clab{icond};
    data_struct.xlabel = 'ISI (ms)';
    data_struct.ylabel = 'statistical frequency';
end
% we can only scale now
if doScale
    maxy = max(H(:));
    for i=1:Ncond
        h = axh(i);
        ylim(h, [0 maxy]);
    end
end
Xlabels(Bh,'ISI (ms)','fontsize',10);
Ylabels(Lh,'statistical frequency','fontsize',10);
% axes(axh(end));
set(gcf,'CurrentAxes',axh(end));
text(0.1, 0.5, IDstring(D, 'full'), 'fontsize', 12, 'fontweight', 'bold','interpreter','none');
if nargout<1, clear H ; end % suppress unwanted echoing
        
function [T,G] = local_ParamGUI
% Returns the GUI for specifying the analysis parameters.
P = GUIpanel('FOISI','');
iCond = ParamQuery('iCond', 'iCond:', '0', '', 'integer',...
    'Condition indices for which to calculate the ISI histogram. 0 means: all conditions.', 20);
Nbin = ParamQuery('Nbin', '# bins:', '50', '', 'posint',...
    'Number of bins used for computing the ISI histogram.', 1);
Anwin = ParamQuery('Anwin', 'analysis window:', 'burstdur', '', 'anwin',...
    'Analysis window (in ms) [t0 t1] re the stimulus onset. The string "burstdur" means [0 t], in which t is the burst duration of the stimulus.');
viewport = ParamQuery('viewport', 'viewport:', 'burstdur', '', 'anwin',...
    'Actual view on the histogram (in ms) [t0 t1]. The string "burstdur" means [0 t], in which t is the burst duration of the stimulus.');
scaling = ParamQuery('scaling','scaling:','',{'auto','constrained'},'',...
    'Click to toggle between histogram scaling options.', 20);
% ParamOrder = ParamQuery('ParamOrder', 'param order:', '', {'[1 2]','[2 1]'}, 'posint',...
%     'Order of independent parameters when sorting [1 2] = "Fastest varied" = "Fastest varied". [2 1] = conversely.', 10);
% SortOrder = ParamQuery('SortOrder', 'sorting order:', '0 0', '', 'integer',...
%     'Sorting order of corresponding independent parameters. (-1,0,1)=(descending, as visited, ascending)',10);
P = add(P, iCond);
P = add(P, Nbin, below(iCond));
P = add(P, Anwin, below(Nbin));
P = add(P, viewport, below(Anwin));
P = add(P, scaling, below(viewport));
% P = add(P, ParamOrder, below(Anwin));
% P = add(P, SortOrder, below(ParamOrder));
P = marginalize(P,[4 4]);
G = GUIpiece([mfilename '_parameters'],[],[0 0],[10 10]);
G = add(G,P);
G = marginalize(G,[10 10]);
% list all parameters in a struct
T = VoidStruct('iCond/Nbin/Anwin/viewport/scaling');
