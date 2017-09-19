function T = IOfun(DS, Chan, iX, varargin);
% dataset/IOfun - input/output curves from RF data
%    IOfun(DS, ADchan) computes and displays I/O curves for each of
%    the frequencies used. ADchan is the channel spec of the analog channel.
%    See dataset/anachan. iCond is the stimulus condition. Default iCond 
%    is 0, meaning all conditions.
%
%    IOfun(DS, ADchan, iX) selects only the x-values (frequencies for RF)
%    with indices iX. Default is iX=0, meaning all frequencies used.
%
%    IOfun(DS, Chan, iX, 'prop1', val1, ...) uses non-default
%    analysisparameters. See Magn_phase_plot for a description.
%
%    The underlying analysis is delegated to dataset/magn_phase_plot.
%  
%    See also dataset/powerspec, dataset/magn_phase_plot.

[iX] = arginDefaults('iX', 0); % default: all freqs
if numel(DS)>1, % must all be sngle-cond datasets; combine them
    error('Multiple DS not allowed yet.');
    % look at dataset/apple how to implement
end

% ========single DS from here================
if ~isequal('RF', upper(stimtype(DS))),
    error('Dataset is not a RC dataset.');
end

% delegated the real work to magn_phase_plot
S = Magn_Phase_Plot(DS, Chan, varargin{:});
% cast Freq, Magn_dB and Phase in X(iSPL, iFreq) matrix form
Freq = [S.Freq].'; 
Magn_dB = [S.Magn_dB].'; 
Phase_cycle = [S.Phase_cycle].';
% also put SPL in this matrix form
SPL = [];
for iSPL = 1:numel(S),
    SPL = [SPL; S(iSPL).Yval.'];
end
% apply iX selection 
if ~isequal(0, iX),
    [SPL Freq Magn_dB Phase_cycle] = deal(SPL(:,iX), Freq(:,iX), Magn_dB(:,iX), Phase_cycle(:,iX));
end
X = S(1).IndepVarX; Y = S(1).IndepVarY; % note that we will swap X&Y re Magn_phase_plot output
Xlab = [Y.ParName ' (' Y.ParUnit ')'];
[dum, Un] = conversionfactor(anachan(DS, Chan));
Ylab1 = ['Magn (dB re 1 ' Un ')'];
Ylab2 = 'Phase (Cycle)';
TitleString = IDstring(DS, 'full');
LegendStr = UnitedLegend(round(Freq(1,:)), X.FormatString);
T1 = structpart(S(1), {'ExpName'    'iDataset'    'Dataset_ID'    'Param'    'IndepVarX'    'IndepVarY'});
T2 = CollectInStruct(SPL, Freq, Magn_dB, Phase_cycle);
T3 = CollectInStruct(TitleString, Xlab, Ylab1, Ylab2, LegendStr);
H = S(1).H; % parameter description
T = structJoin(T1, '-mag_phase', T2, '-plot_accessories', T3, '-param_descr', H);
if nargout<1, 
    local_plot(T);
end


%============================================
function ha = local_plot(T);
set(gcf,'units', 'normalized', 'position', [0.333 0.18 0.534 0.722]);
if nargin<2, clf; CLR = lico(1); end
ha(1) = subplot(2,1,1); 
plot(T.SPL, T.Magn_dB, 'o-'); 
set(ha(1), 'xtick', -300:10:300, 'ytick', -300:10:300, 'fontsize', 12);
grid on;
ylabel(T.Ylab1, 'fontsize', 14);
legend(T.LegendStr, 'location', 'northwest', 'fontsize', 8); 
title(T.TitleString);
ha(2) = subplot(2,1,2); 
plot(T.SPL, cunwrap(T.Phase_cycle), 'o-'); 
set(ha(2), 'xtick', -300:10:300, 'ytick', -100:0.5:100, 'fontsize', 12);
xlabel(T.Xlab, 'fontsize', 14);
ylabel(T.Ylab2, 'fontsize', 14);
grid on;





