function [XFit, YFit, Env, Param] = gaborfit(X, Y, NewRange, SampleRate, varargin)

%B. Van de Sande 25-03-2003

%% ---------------- CHANGELOG -----------------------
%  Tue Apr 24 2012  Abel   
%   - re-write to support initial estimates for fitting
%[Tue Jul 14 2015 (Abel)]: 
% - updated syntax for optimset ('Algorithm','levenberg-marquardt')


%% ---------------- Default parameters ---------------
defParam = struct();
%Constants for Gabor function
% [ Amplitude, 
%   Position of enveloppe-maximum, ...
%   Width of enveloppe,...
%   Shape of enveloppe (small values give accentuation of the central peak, large values give a box-like shape),... 
%   Frequency,...
%   Phase-shift ]
% format: [ LowerBound Estimate UpperBound ]
% original values (set by Bram):
% 	CEstimate   = [100,     0,  750,   1/1000, 0];
% 	LowerBounds = [25,  -1000,  100, 0.5/1000, -pi];
% 	UpperBounds = [200, +1000, 5000,  50/1000, +pi];

defParam.ampl = [25 100 200];
defParam.max = [-1000 0 1000];
defParam.width = [100 750 5000];
defParam.freq = [0.5/1000 1/1000 50/1000];
%GET pH min/max in order
defParam.ph = [-pi 0 pi];

%% ---------------- Main function --------------------
if nargin > 4
	param = getarguments(defParam, varargin);
	
	%reset NaN's to default values
	fields = fieldnames(defParam);
	for n = 1:length(fields)
		if isnan(param.(fields{n}))
			param.(fields{n}) = defParam.(fields{n});
		end
	end
	
	%populate upper and lower bounds if needed 
	% - amplitude, gabor function width, gabor function freq
	fields = {'ampl','width','freq'};
	for n = 1:length(fields)
		if length(param.(fields{n})) ~= 3
			lowFact = safeFact_(defParam.(fields{n})(1), defParam.(fields{n})(2));
			highFact = safeFact_(defParam.(fields{n})(3), defParam.(fields{n})(2));
			param.(fields{n}) = [ param.(fields{n})*lowFact, param.(fields{n}), param.(fields{n})*highFact ];
		end
	end
	% - position of env max
	if length(param.max) ~=3
		param.max = [ param.max - (defParam.max(2) - defParam.max(1)),...
			param.max,...
			param.max + (defParam.max(2) + defParam.max(3)) ];
	end
	% - phase of gabor function
	if length(param.ph) ~=3
		param.ph = [ param.ph - (2*pi),...
			param.ph,...
			param.ph + (2*pi) ];
	end		
else
	param = defParam;
end

CEstimate   = [param.ampl(2), param.max(2), param.width(2), param.freq(2), param.ph(2)];
LowerBounds = [param.ampl(1), param.max(1), param.width(1), param.freq(1), param.ph(1)];
UpperBounds = [param.ampl(3), param.max(3), param.width(3), param.freq(3), param.ph(3)];

% 
% disp('--------------------');
% param
% disp('--------------------');
% 
% hold on;
% plot(X,Y, 'g');
    
%Least Square Fit toepassen ...
Options = optimset('Display', 'off','Algorithm','levenberg-marquardt'); %'MaxFunEvals', realmax
% [c, Residual, Dummy, Converged] = lsqcurvefit(@gaborfunc, CEstimate, X, Y, LowerBounds, UpperBounds, Options);
[c, Residual, Dummy, Converged] = lsqcurvefit(@gaborfunc, CEstimate, X, Y, [], [], Options);

if Converged,
    XFit = NewRange(1):(1/SampleRate):NewRange(2);
    YFit = gaborfunc(c, XFit);
    
    %Parameters samenstellen ...
    A  = c(1); %Amplitude ...
    DC = 0;    %DC-value ...
    
    EnvMax   = c(2); %Position of enveloppe-maximum ...
    EnvWidth = c(3); %Width of enveloppe ...
    EnvShape = 2;    %Shape of enveloppe (small values give accentuation of the central peak, large values give a box-like shape)...
    
    Freq = c(4); %Frequency ...
    Ph   = c(5); %Phase-shift ...
    
    %Fractie van variantie die fit in rekening brengt ...
    AccFraction = getaccfrac(@gaborfunc, c, X, Y);
    
    Param = CollectInStruct(A, DC, EnvMax, EnvWidth, EnvShape, Freq, Ph, AccFraction);
    
    %Enveloppe bepalen ...
    Env = A * exp(-((abs(XFit-EnvMax)/EnvWidth).^EnvShape));
else
    NElem = (NewRange(2) - NewRange(1)) * SampleRate;
    [XFit, YFit, Env] = deal(repmat(NaN, 1, NElem));
    [A, DC, EnvMax, EnvWidth, EnvShape, Freq, Ph, AccFraction] = deal(NaN);
    
    Param = CollectInStruct(A, DC, EnvMax, EnvWidth, EnvShape, Freq, Ph, AccFraction);
end

%% ---------------- Local functions ------------------
function divRes = safeFact_(t, n)
if t == 0
	t = 0.01;
end
if n == 0
	n = 0.01;
end
divRes = t/n;


