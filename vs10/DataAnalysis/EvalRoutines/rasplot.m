function Hdl = rasplot(ds, varargin)
%RASPLOT  creates raster plot
%   RASPLOT(ds) creates raster plot of dataset ds. The plot is created in the current axis.
%   If no current axis object exists, one is created.
%   Hdl = RASPLOT(ds, iSub, AxHdl) only plots the requested subsequences and uses the handle
%   of the axis object suppied in AxHdl. Optionally, the handle of the axis object is returned.
%INPUT
%   colors  = {'r', 'g', 'b', 'c', 'm', 'y', 'k'};  %Colors per subsequence
%   isubseqs = 'all';                               %Dataset subsequences to be plotted
%   alpha = 1;                                      %Transparency of plot [0 -> 1] with 0 invisible 
%   axhdl = [];                                     %Optional handle of axes to plot
%                                                    on
%OUTPUT
% Hdl		%Handle of 

%B. Van de Sande 24-05-2004

%% ---------------- CHANGELOG ------------------------
%  Thu Jul 7 2011  Abel   
%  - Added support for plotting transparent RAS plots
%  - Added return default options functionality 
%% ---------------- Default parameters ---------------
Defaults.colors  = {'r', 'g', 'b', 'c', 'm', 'y', 'k'};
Defaults.isubseqs = 'all';
Defaults.alpha = 1;
Defaults.axhdl = [];

%% ---------------- Main function --------------------
% Return if no arguments
if nargin < 1
	Hdl = Defaults;
	return;
end 

% Get params
param = getarguments(Defaults, varargin);
Colors = param.colors;
if ~isa(ds, 'dataset')
	printhelp(Defaults, 'First argument should be a dataset');
end

if strcmpi(param.isubseqs, 'all')
	iSubSeqs = 1:ds.Stim.Presentation.Nrep;
else
	iSubSeqs = param.isubseqs;
end

if isempty(param.axhdl)
	AxHdl    = gca;
else 
	AxHdl = param.axhdl;
end

% Set colors, repeat single color if needed
NColors = length(Colors);
if NColors == 1
	Colors = repmat(Colors, length(iSubSeqs));
	NColors = length(iSubSeqs);
end

%Plotting raster ...
SptCell = spiketimes(ds);
NRep    = ds.Stim.Presentation.Nrep;
ColIdx  = 0; %Color index ...
Hdl     = [];
for nSubSeq = iSubSeqs(:)'
    X = [];
    Y = [];
    ColIdx = mod(ColIdx, NColors) + 1;
    yy = linspace(nSubSeq-0.4, nSubSeq+0.4, NRep+1);
    for n = 1:NRep
        Spks = SptCell{nSubSeq, n}; 
        NSpk = length(Spks);
        if (NSpk > 0)
            X  = [X VectorZip(Spks, Spks , Spks+NaN)];
            y1 = yy(n)+0*Spks;   %same size as x
            y2 = yy(n+1)+0*Spks; %same size as x
            Y  = [Y, VectorZip(y1, y2, y2)];
        end
	end
	%Plot using patch(), needed for transparency 
	LnHdl = patch(X, Y, Colors{ColIdx}, 'EdgeColor', Colors{ColIdx}, 'EdgeAlpha', param.alpha);
    if ~isempty(LnHdl)
        Hdl = [Hdl, LnHdl];
    end
end

%Setting axis properties ...
%Indepval vector is invalid for incomplete BFS-dataset collected before August 2002 ...
Indepval = ds.Stim.(ds.Stim.Presentation.X.FieldName);
set(AxHdl, 'YTick', iSubSeqs, 'YTicklabel', eval(['{' num2sstr(Indepval(iSubSeqs))  '}']), ...
    'ylim', [min(iSubSeqs)-0.5 max(iSubSeqs)+0.5]);
title('Raster Plot', 'fontsize', 12);
xlabel('Time (ms)');
ylabel(sprintf('%s (%s)', ds.Stim.Presentation.X.ParName, ds.Stim.Presentation.X.ParUnit));

%% ---------------- Local functions ------------------

