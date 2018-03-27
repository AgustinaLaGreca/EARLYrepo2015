function ah=plot(EC, varargin);
% Earcalib/plot - plot Earcalib object
%   plot(PTC) plots Earcalib object PTC, displaying the Magnitude and
%   phase transfers of the Probe & Cavity recordings in one pair of graphs
%   and those of the  Probe-to-Cavity transfer in aonther pair. Plot
%   returns graphics handles to the four axes systems.
%
%   plot(ah, PTC) plots to axes with handles ah(1) to ah(4).
%
%   See also Transfer/Plot.

figh = gcf;
figh.IntegerHandle = 'on';
figh = figh.Number;
set(findobj(figh,'type', 'line'), 'linewidth',1);

if ~isa(EC, 'earcalib'), % plot(ah, T, ...) syntax. Rearrange input args
    if nargin<2 || ~isa(varargin{1}, 'earcalib'),
        error('Either the first or the second input arg of Earcalib/Plot must be Earcalib object.');
    end
    if any(~istypedhandle(EC,'axes')) || ~isequal(4,numel(EC)),
        error('First input arg of Earcalib/Plot must be either a Earcalib object or a quadruple of axes handles.');
    end
    % now rearrange the input args
    ah = EC;
    EC = varargin{1};
    varargin = varargin(2:end);
else, % plot(T, ...) syntax. Provide own axes handles.
    set(figh,'units', 'normalized', 'position', [0.324 0.198 0.653 0.716]);
    ah = [subplot(2,2,1) subplot(2,2,3) subplot(2,2,2) subplot(2,2,4)];
end
ah = double(ah);
set(ah,'visib','off'); % plotting will make them reappear
Trf1 =EC.Transfer(1); Trf2 =EC.Transfer(2);

icol1 = 1+getGUIdata(figh,'ColorIndex_1',0);
LegStr = [getGUIdata(figh,'LegStr',{}) [upper(name(EC.Experiment)) ' calib # ' num2str(EC.iCalib)]];
setGUIdata(figh,'LegStr',cellify(LegStr));
if ~isvoid(Trf1),
    setGUIdata(figh,'ColorIndex_1',icol1);
    mFactor = unique(MaxAbsDA(Trf1))/sqrt(2); % factor from 1 Volt RMS to a peak Voltage of maxabsDA(Trf1)
    set(ah(1:2),'visib','on')
    plot(ah(1:2), mFactor*Trf1, lico(icol1), 'linewidth', 1.5, varargin{:});
    axes(ah(1)); xlog125(fmin(Trf1)/1e3, fmax(Trf1)/1e3);
    set(ah(2), 'xlim', [0, fmax(Trf1)]/1e3);
    ylabel(ah(1),'Max Tone Level (dB SPL)');
    legend(ah(1), LegStr, 'location', 'southwest');
end

icol2 = 1+getGUIdata(figh,'ColorIndex_2',0);
if ~isvoid(Trf2),
    setGUIdata(figh,'ColorIndex_2',icol2);
    mFactor = unique(MaxAbsDA(Trf2))/sqrt(2); % factor from 1 Volt RMS to a peak Voltage of maxabsDA(Trf1)
    set(ah(3:4),'visib','on');
    plot(ah(3:4), mFactor*Trf2, lico(icol2), 'linewidth', 1.5, varargin{:});
    axes(ah(3)); xlog125(fmin(Trf2)/1e3, fmax(Trf2)/1e3);
    set(ah(4), 'xlim', [0, fmax(Trf2)]/1e3);
    ylabel(ah(3),'Max Tone Level (dB SPL)');
    legend(ah(3), LegStr, 'location', 'southwest');
end
if nargout<1, clear ah; end

