clear
%last number is the index in the index of the stimulus in the experiment
D = read(dataset,'F14535',3);

%47
%55
%65,67

start_test=500;
end_test = 700;



stimparapm = D.stimparam;

iCond=stimparapm.Presentation.iCond;
iRep = stimparapm.Presentation.iRep;

Xname = stimparapm.Presentation.X.FieldName;
Xval = stimparapm.Presentation.X.PlotVal;
nX=length(unique(Xval));

Yname = stimparapm.Presentation.Y.FieldName;
Yval = stimparapm.Presentation.Y.PlotVal;
nY=length(unique(Yval));

% NX,NY = stimparapm.Ncond_XY;
Nconds = prod(stimparapm.Ncond_XY);
spks = D.spiketimes;
Nreps = stimparapm.Nrep

count_mat = zeros(nY,nX);

bins = 0:10:2000;
lPSTH = length(bins);

count=zeros(Nconds);
PSTHs = zeros(Nconds,lPSTH);
for icond=1:Nconds
    for irep=1:Nreps
        train = spks{icond,irep};
%         PSTHs(icond,:)=PSTHs(icond,:)+hist(train,bins);
        ispk = intersect(find(train>=start_test),find(train<=end_test));
        train = train(ispk);
        PSTHs(icond,:)=PSTHs(icond,:)+hist(train,bins);
        count(icond) = count(icond)+length(train);
    end
end

icond=1
for iY=1:nY
    for iX=1:nX
        count_mat(iY,iX)=count(icond);
        icond=icond+1;
    end
end


% figure()

figure()
subplot(211)
plot(PSTHs')
subplot(212)
plot(Xval,count_mat)

% imagesc(count_mat)

