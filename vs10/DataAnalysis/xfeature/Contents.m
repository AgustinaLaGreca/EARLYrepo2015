% XFEATURE
%
% TK AFB data analysis
%   deristats             - view statistics of time-derivative of recordings
%   Xfeature_sandbox      - test ground for Xfeature functionality derived from ...
%
% Chunks and snippets.
%   xchunks               - basic chunk extraction from recording trace
%   EvalChunks            - evaluate chunks of a signal with respect to given property
%   getchunk              - get single segment from signal
%   ChunkClick            - list chunk that is clicked 
%   BasicChunkPlot        - plot samples and highlight chunks
%   chunkPlot             - plot or highlight chunks
%   getsnips              - get snippets from recorded TK ABF waveform
%   categorize            - divide array elements over bins
%
% Plots.
%   deriplot              - plot deristat output
%   chunkPlot             - plot or highlight chunks
%   CatPlot               - plot categorized columns of matrix and derived templates
%   sscatter              - slow scatter plot
%   qscatter              - quick scatter plot
%   XcycleHist            - cycle histogram based on AP/EPSP output of deristats
%   TracePlotInterface    - enable keyboard navigation for plots of recording traces
%   BasicChunkPlot        - plot samples and highlight chunks
%
% Helpers.
%   catsam                - concatenate all samples of TK ABF file
%   smoothen              - smooth data by convolution with a hamming window
%   XgetStimprops         - stimulus properties in a struct
%   localmax              - determine local maxima of array
%   mostRecent            - for given instant, determine most recent instant from set
%   vectorstrength        - vector strength of train of event times
%   MultiMedian           - medians across various columns of a matrix.
%   NormDist              - nonsymmetric normalized distance
%   locateSegment         - locate segment within array
%   testPCAwave           - test waveforms for PCA analysis
%   derivative            - estimated time-derivative of sampled signal
%
% Old stuff.
%   xspikeAlign           - time-align individual events with their template
%   xspikeplot            - plot spikes extracted by xspikeX.
%   xspikes               - extract spikes from recorded traces
%   templot               - aligned plot of action potential templates
%   SpikeMetrics          - evaluate spike metrics: max EPSP rate, AP peak, etc.
%   PlotSpikeMetrics      - plot trace annotated with spike metrics 
%   tchebi                - event extraction using tchebishev polynomials
%   SPMcycleHist          - cycle histogram from spike metrics output
%   Mtemplate             - event templates based on spike metrics
%   getevent              - get single event from xspike output
%   getABFfromSpikeStruct - retrieve ABF data that produced xspikes output S
%   old_deristats         - view statistics of time-derivative of recordings 