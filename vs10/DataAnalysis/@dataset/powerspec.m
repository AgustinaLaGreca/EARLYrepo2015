function PP = powerspec(DS, ADchan, iCond, iRep, Twin, varargin);
% dataset/powerspec - power spectrum of analog recording
%    powerspec(DS, ADchan, iCond, iRep, Twin) plots the power spectrum of the
%    analog responses in analog channel ADchan of dataset D. 
%    For specifying ADchan, see dataset/anachan. Twin is the analysis 
%    window applied to the recorded waveforms. The default Twin is 
%    "burstdur". If multiple repetitions are present for the specified 
%    condition, their average is taken (see dataset/anamean). One or several
%    individual repetitions may be selected bny iRep; the default iRep is 0, 
%    meaning all reps.
%
%    powerspec(DS, ADchan, iCond, iRep, Twin, ...) passes any trailing input
%    arguments to PLOT.
%
%    P = powerspec(D, ADchan, Twin, Cdelay) suppresses plotting and returns
%    the power spectrum in struct P.
%
%    See also dataset/Magn_Phase_Plot, dataset/anamean, dataset/anachan.


if nargin<4, [iRep, Twin] = arginDefaults('iRep/Twin', 0, []); % all reps / burst only
elseif nargin<5, Twin = []; end;
if isempty(Twin), Twin = 'BurstDur'; end

if ischar(Twin), % if Twin resembles 'BurstDur', turn it to 'BurstDur' exactly. Same for 'RepDur'
    Twin = keywordMatch(Twin, {'BurstDur' 'RepDur'});
end


if numel(iCond)>1, % recursive call
    for ii=1:numel(iCond),
        if nargout>1,
            PP(ii) = powerspec(DS, ADchan, iCond(ii), iRep, Twin);
        else,
            powerspec(DS, ADchan, iCond(ii), iRep, Twin, lico(ii));
            legend(gca, 'off');
        end
    end
    if nargout<1, % i.e. when plotting
        H = findobj(gcf, 'type', 'line', 'marker', 'none');
        H = flipud(H(:));
        legend(H, cellify(CondLabel(DS, iCond)), 'location', 'northwest'); 
    end
    return;
end

%============single stim condition from here==============

IDstr = IDstring(DS, 'full');
ExpName = expname(DS);
Irec = DS.ID.iDataset;
P = CollectInStruct(IDstr, ExpName, Irec, ADchan, iCond, iRep, Twin, '-');

if isequal('BurstDur', Twin), 
    Twin = [0 max(burstdur(DS, iCond))];
elseif isequal('RepDur', Twin), 
    Twin = [0 repdur(DS, iCond)];
end

% get analog recording
[D, dt, t0, DataType] = anamean(DS, ADchan, iCond, iRep, Twin);
Nsam = numel(D);
[ScalFac, Unit] = conversionfactor(anachan(DS, ADchan));
D = D*ScalFac; % rec in true units (not volts of the ADC)


P.df = Fsam(anachan(DS, ADchan))/Nsam; % sample freq in Hz
P.Power_db = A2dB(2*abs(fft(D.*hann(Nsam)/Nsam))) - P2dB(P.df); % spectrum level
NsamPos = floor(Nsam/2); % max index of pos freqs components
P.Power_db = P.Power_db(1:NsamPos);
if isequal('Pa', Unit),
    P.Power_db = P.Power_db - A2dB(20e-6); % 1 Pa = "94 dB"
    P.Yunit = 'dB SPL/Hz';
else,
    P.Yunit = ['dB re 1 ' Unit];
end
P.df_unit = 'Hz';

% is there a "Fcar" in the stim def of DS? If so, know how to mark it.
[P.Fstim, idx] = local_stimfreq(DS, iCond, P.df);
P.Power_dB_Fstim = P.Power_db(idx);
if nargout<1, % plot
    set(gcf,'units', 'normalized', 'position', [0.333 0.301 0.58 0.601]);
    if numel(varargin)>0, % add to existing plot
        hl = xdplot(P.df/1e3, P.Power_db, varargin{:});
    else,
        hl = xdplot(P.df/1e3, P.Power_db);
    end
    xplot(P.Fstim/1e3, P.Power_dB_Fstim, varargin{:}, 'marker', '*','color','k');
    xplot(P.Fstim/1e3, P.Power_dB_Fstim, varargin{:}, 'marker', 'o','color','k');
    xlog125;
    xlabel('frequency (kHz)','fontsize',10);
    ylabel(['power (' P.Yunit ')'],'fontsize',10);
    grid on
    legend(hl, CondLabel(DS, iCond));
    title([IDstring(DS, 'full') '  AD-'  num2str(ADchan) ' (' dataType(anachan(DS,ADchan)) ')'], 'fontsize', 12, 'fontweight', 'bold', 'interpreter', 'none');
else, % return P
    PP = P;
end

function [Fstim, idx] = local_stimfreq(DS, iCond, df);
[Fstim, idx] = deal(nan,1);
Stim = DS.Stim;
if isfield(Stim, 'Fcar'),
    Fstim = Stim.Fcar; %  in Hz
    if size(Fstim,1)>1, % Fcar changes w condition
        Fstim = Fstim(iCond, :);
    end
    if ~isnan(Fstim), idx = 1+round(Fstim/df); end
end











