function [phase_values,freq]=plotPhase(W, varargin);
% Waveform/plotPhase - plot phase of frequency spectrum of waveform object
%     plotPhase(W, ...) plots phase spectrum of waveform object

% JV 9/7/2015

if size(W,1)>1, error('Cannot plot multiple waveforms unless L/R pair'); end

%set(gcf,'units', 'normalized', 'position', [0.6225 0.0775 0.35 0.35])

CLR = get(0,'defaultAxesColorOrder'); NCL = size(CLR,1);

Nchan = size(W,2);
Fsam = W(1).Fsam; % sample frequency in s
if isequal('replace', get(gca,'NextPlot')), cla; end
LegStr = {};
if ~isempty(varargin) 
    if strcmpi(varargin{1},'diff')
        if strcmpi(W(1).DAchan,'r')
            ichan_r = 1;
            ichan_l = 2;
        else
            ichan_r = 1;
            ichan_l = 2;
        end
        x_left = samples(W(ichan_l));
        x_right = samples(W(ichan_r));
        % single sided amplitude spectrum
        N = 2^nextpow2(size(x_right-x_left,1));
        df = Fsam/N;
        X = angle(fft(x_right-x_left,N))/(2*pi);
        if nargout > 0
            phase_values = X(1:N/2+1);
            freq = df;
        end
        xdplot(df,X(1:N/2+1),'color', 'r');
        LegStr{end+1} = 'Diff (R - L)';    
    end
else       
    for ichan=1:Nchan,
        x = samples(W(ichan));
        switch W(ichan).DAchan,
            case 'L', clr = CLR(1,:);
            case 'R', clr = CLR(2,:);
            otherwise, clr = CLR(1+rem(ichan-1,NCL),:);
        end
        % single sided amplitude spectrum
        N = 2^nextpow2(size(x,1));
        df = Fsam/N;
        X = fft(x,N);
        threshold = max(X)/10;
        X(abs(X)<threshold) = 0; 
        X = angle(X)/(2*pi);
        if nargout > 0
            phase_values = X(1:N/2+1);
            freq = df;
        end
        xdplot(df,X(1:N/2+1),'color', clr, varargin{:});
        LegStr{end+1} = W(ichan).DAchan;    
    end
end
xlim('auto');
XL = xlim;
%xlim([mean(XL)+1.1*(XL-mean(XL))]);
%   Fsam: 3.2552e+004
%        DAchan: 'R'
%     MaxMagSam: 3.1623e-005
%         Param: [1x1 struct]
%       Samples: {[0]  [127x1 double]  [96x1 double]  [65x1 double]  [0]}
%          Nrep: [1 50 1 1 6509]
%

if ~isempty(varargin) 
    if ~strcmpi(varargin{1},'diff')
        xlabel('Frequency (Hz)');
        ylabel('Phase (cycles)');
    end
else
    xlabel('Frequency (Hz)');
    ylabel('Phase (cycles)');
    legend(LegStr{:});
end
