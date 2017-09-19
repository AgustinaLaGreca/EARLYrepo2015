function cHist = SGSR_cyclehist(varargin)
%CYCLEHIST	calculate cyclehistogram.
%   cHist = CYCLEHIST(X, Y, AnDur, NReps, BinFreq, {NBin, CacheParam, CacheFileName})
%   Input parameters:
%   X       : vector with time-axis in milliseconds
%   Y       : vector with spike rate in spk/sec
%   AnDur   : duration of analysis window in ms
%   NReps   : repetitions use in calculating the correlogram
%   BinFreq : frequency at which to calculate the cyclehistogram
%   NBin    : number of bins in the cyclehistogram {64}
%   Attention! Calculating these cyclehistograms is time-consuming, these calculations can therefore be cached.
%   If supplied as argument, the results are saved in the file CacheFileName. If no directory is given the current
%   directory is assumed. CacheParam is a structure with uniquely identifies the calculated data.
%   Output parameters:
%   cHist is a structure with following fields:
%             cHist.X       : from 0 to 1 in NBin steps
%             cHist.Y       : coincidence probability
%             cHist.R       : vector strength magnitude
%             cHist.Ph      : vector strength phase
%             cHist.pRaySig : Rayleigh significance
%             cHist.N       : number of spikes
%
%	cHist = CYCLEHIST(ds, NSub, {BinFreq, NBin, AOn, AOff})
%   Input Parameters:
%   ds      : dataset-object
%   NSub    : subsequence number
%   BinFreq : frequency at which to calculate the cyclehistogram {-1 = Stimulus Frequency}
%   NBin    : number of bins in the cyclehistogram {64}
%   AOn     : onset of analysiswindow {0ms}
%   AOff    : end of analysiswindow {-1 = Stimulus Duration}

%B. Van de Sande 15-12-2004

%% ---------------- CHANGELOG -----------------------
%  Tue Jul 12 2011  Abel   
%   - Removed bug for matlab2009 and above 

if isempty(varargin), error('Wrong number of input arguments.'); end

if isa(varargin{1}, 'double')
	switch nargin
		case 5
			X             = varargin{1};
			Y             = varargin{2};
			AnDur         = varargin{3};
			NReps         = varargin{4};
			BinFreq       = varargin{5};
			NBin          = 64;
			CacheParam    = struct([]);
			CacheFileName = '';
		case 6
			X             = varargin{1};
			Y             = varargin{2};
			AnDur         = varargin{3};
			NReps         = varargin{4};
			BinFreq       = varargin{5};
			NBin          = varargin{6};
			CacheParam    = struct([]);
			CacheFileName = '';
		case 8
			X             = varargin{1};
			Y             = varargin{2};
			AnDur         = varargin{3};
			NReps         = varargin{4};
			BinFreq       = varargin{5};
			NBin          = varargin{6};
			CacheParam    = varargin{7};
			CacheFileName = varargin{8};
		otherwise, error('Wrong number of input arguments.'); end
		
		if rem(NBin, 2), error('NBin must be even number.'); end
		if isempty(Y) | all(Y == 0), cHist = struct('X', linspace(0,1, NBin), 'Y', repmat(NaN, 1, NBin), 'R', NaN, 'Ph', NaN, 'pRaySig', NaN, 'N', NaN); return; end
		
		if ~isempty(CacheFileName),
			[Path, FileName, FileExt] = fileparts(CacheFileName);
			if isempty(Path),    Path = cwd; end
			if isempty(FileExt), FileExt = ['.' mfilename '.cache']; end
			if isempty(FileName), error('Invalid cache filename.'); end
			CacheFileName = fullfile(Path, [FileName FileExt]);
		end
		if ~isempty(CacheParam), CacheParam = lowerFields(rofields(structcat(CacheParam, CollectInStruct(AnDur, NReps, BinFreq, NBin)), 'alpha')); end
		
		Data = local_FromCache(CacheFileName, CacheParam);
		if ~isempty(Data), unpackstruct(Data);
		else
			%Basisgegevens verzamelen ...
			Period   = 1000/BinFreq;   %Duur van periode in ms
			BinWidth = X(2) - X(1);    %Duur van een bin in ms
			MaxLag   = X(end);         %Maximum tijdsverschil in ms
			
			%Nagaan hoeveel cycli in tijdsinterval gaan, indien laatste cyclus niet volledig in tijdsinterval past
			%dan dit laten wegvallen ...
			NPeriodML = floor((MaxLag + (BinWidth/2)) / Period);
			NewMaxLag = Period * NPeriodML;
			NPeriod   = NPeriodML * 2;
			NBinPeriod = Period/BinWidth;
			PeriodRatio = NBinPeriod/NBin;
			
			%Denormaliseren ...
			Y = round(Y *NReps *AnDur /1000);
			%Unbinnen: de curve beschouwen als een histogram en dit omzetten naar een spiketrain...
			Spt = []; 
            for ind = unique(Y), 
                indices = find(Y == ind);
                if ~isempty(indices)  
                	Spt = [Spt repmat(X(find(Y == ind)), [1, ind])]; 
                end
            end
            %Dittering zodat de spikes toch verspreid liggen en dus niet samengeklusterd zijn thv van bincenter ...
			if ~isempty(Spt)
                Spt = Spt + ((rand(1,length(Spt))-0.5) * BinWidth);
            end
            %Spiketimes die buiten de grenzen van een het nieuwe tijdsinterval vallen weglaten, deze maken immers slechts deel
			%uit van een onvolledige cyclus ...
			Spt = Spt(find((Spt >= -NewMaxLag) & (Spt < NewMaxLag)));
			%De negatieve spiketimes alle positief maken ... Nul wordt automatisch referentie ...
			Spt = [Spt(find(Spt >= 0)) (Spt(find(Spt < 0)) + NewMaxLag)];
			
			%Herbinnen volgens de periode van de stimulus ...
			Cycle = 0:1/NBin:1;
			NCo = histc( rem(Spt, Period)/Period, Cycle);
            if ~isempty(NCo)
                NCo(end) = [];
            end
            
			N = length(Spt);
			Data = CollectInStruct(NCo, NPeriod, PeriodRatio, N);
			local_ToCache(CacheFileName, CacheParam, Data);
		end
		
		X = 1/(2*NBin):1/NBin:(1-1/(2*NBin));
		%Normaliseren naar spikes/sec ...
		Y = 1000 * NCo /PeriodRatio /AnDur /NPeriod /NReps;
        if isempty(Y)
            x_dim = size(Y,1)-1;
            y_dim = size(Y,2)+1;
            Y=zeros(y_dim,x_dim);
        end
		%Vector strength berekenen
		% 1) Omzetten van polaire naar cartesiaanse coordinaten
		PhRad = X * (2*pi);
		x     = Y .* cos(PhRad);
		y     = Y .* sin(PhRad);
		% 2) Vectorsom nemen
		x = sum(x);
		y = sum(y);
		% 3)Grootte berekenen via stelling van Pytagoras en
		%	 normaliseren
		R = sqrt(x^2+y^2) / sum(Y);
		% 4) Phase van resulterende vector berekenen door terug te keren
		%    naar polaire coordinaten
		Ph = atan2(y,x) / (2*pi);
		if Ph < 0, Ph = Ph +1; end
		
		%Rayleigh significantie berekenen
		pRaySig = RayleighSign(R, N);
		
		if nargout > 0, cHist = CollectInStruct(X, Y, R, Ph, pRaySig, N); end
elseif isa(varargin{1}, 'dataset')
	switch length(varargin)
		case 2
			ds       = varargin{1};
			NSub     = varargin{2};
			BinFreq  = -1;
			NBin     = 64;
			AOn      = 0;
			AOff     = -1;
		case 6
			ds       = varargin{1};
			NSub     = varargin{2};
			BinFreq  = varargin{3};
			NBin     = varargin{4};
			AOn      = varargin{5};
			AOff     = varargin{6};
        otherwise
            error('Wrong number of input arguments.'); end
		
		if strncmpi(ds.StimType, 'THR', 3), error('Dataset doesn''t contain spiketime-information.'); end
		if NSub > ds.Stim.Presentation.Ncond , error('Subsequence doesn''t exist or wasn''t recorded.'); end
		
		%Changed on 19-11-2003 ...
		if AOff == -1, AOff = max(ds.Stim.BurstDur); end
        if BinFreq < 0
            if strcmp(ds.Stim.DAC,'Both') %Binaural data ...
                if ~all(isnan(ds.Stim.ModFreq)) & ~all(ds.Stim.ModFreq == 0),
                    %Return the column that doesn't contains only NaN's ...
                    idx = find(~all(isnan(ds.Stim.ModFreq), 1));
                    if ~isempty(idx), BinFreq = abs(BinFreq) * ds.Stim.ModFreq(:, idx);
                    else, error('No standard binning frequency for dataset.'); end
                elseif ~all(isnan(ds.Stim.Fcar)) & ~all(ds.Stim.Fcar == 0),
                    %Return the column that doesn't contains only NaN's ...
                    idx = find(~all(isnan(ds.Stim.Fcar), 1));
                    if ~isempty(idx), BinFreq = abs(BinFreq) * ds.Stim.Fcar(:, idx);
                    else, error('No standard binning frequency for dataset.'); end
                else, error('No standard binning frequency for dataset'); 
                end
            else, %Monaural data ...
                %by Abel: Update for matlab >=2009:
                %         ds.Special.ModFreq is not equal to 0 if its
                %         undefined, but to the smallest number known by matlab
                %         (1.4694e-039). We therefore need to include an extra
                %         check in the if-function: ds.Special.ModFreq < 0.001
                if ~all(isnan(ds.Stim.ModFreq)) && ~all(ds.Stim.ModFreq < 0.001)
                    BinFreq = abs(BinFreq) * ds.Stim.ModFreq;
                    % 			elseif ~all(isnan(ds.Special.ModFreq)) && ~all(ds.Special.ModFreq == 0)
                    % 				BinFreq = abs(BinFreq) * ds.Special.ModFreq;
                elseif ~all(isnan(ds.Stim.Fcar)) && ~all(ds.Stim.Fcar < 0.001)
                    BinFreq = abs(BinFreq) * ds.Stim.Fcar;
                    % 			elseif ~all(isnan(ds.Special.CarFreq)) && ~all(ds.Special.CarFreq == 0)
                    % 				BinFreq = abs(BinFreq) * ds.Special.CarFreq;
                else
                    error('No standard binning frequency for dataset');
                end
            end
        end
        %The binning frequency is always a scalar ...
        if (length(BinFreq) ~= 1), BinFreq = BinFreq(NSub); end
		
        Period  = 1000/BinFreq;      %duur van periode in ms
		NPeriod = (AOff-AOn)/Period; %aantal perioden in analysewindow
		
        spt = spiketimes(ds);
		spt = spt{NSub};
		sptime = anwin(cat(2, spt(:)), [AOn AOff]);
		N = length(sptime);
		
		%Indien deze subsequentie geen spiketimes bevat ...
		if isempty(sptime), cHist = struct('X', 1/(2*NBin):1/NBin:(1-1/(2*NBin)), 'Y',zeros(1, NBin), 'R', NaN, 'Ph', NaN, 'pRaySig', NaN, 'N', NaN); return; end
		
		Cycle = 0:1/NBin:1;
		nspt = histc( rem(sptime, Period)/Period, Cycle);
		nspt(end) = [];
		
		X = 1/(2*NBin):1/NBin:(1-1/(2*NBin));
		Y = 1000 * NBin *nspt/(AOff-AOn) /ds.Stim.Presentation.Nrep;
		
		%Vector strength berekenen
		% 1) Omzetten van polaire naar cartesiaanse coordinaten
        if size(X') == size(nspt)
            X=X';
        end
		PhRad = X*2*pi;
		x = nspt .* cos(PhRad);
		y = nspt .* sin(PhRad);
		% 2) Vectorsom nemen
		x = sum(x);
		y = sum(y);
		% 3)Grootte berekenen via stelling van Pytagoras en
		%	 normaliseren
		R = sqrt(x^2+y^2) / length(sptime);
		% 4) Phase van resulterende vector berekenen door terug te keren
		%    naar polaire coordinaten
		Ph = atan2(y,x)/ (2*pi);
		if Ph < 0, Ph = Ph + 1; end
		
		%Rayleigh significantie berekenen
		pRaySig = RayleighSign(R, length(sptime));
		
		if nargout > 0, cHist = CollectInStruct(X, Y, R, Ph, pRaySig, N); end
else, error('Wrong input arguments.'); end

%--------------------------locals------------------------------------
function local_ToCache(FileName, SParam, Data)

if isempty(FileName), return; end

Data.NCo = ZipNumericData(Data.NCo);
PutInHashFile(FileName, SParam, Data, +1021);

function Data = local_FromCache(FileName, SParam)

Data = [];

if isempty(FileName), return; end

Data = GetFromHashFile(FileName, SParam);
if isempty(Data), return; end

Data.NCo = UnZipNumericData(Data.NCo);