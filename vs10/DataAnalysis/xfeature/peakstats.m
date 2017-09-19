function peakstats(dt, X, DT, Nbin);
% peakstats - display stats on waveform peaks
%    peakstats(dt, X, DT, Nbin) displays the distribution of the running maxima
%    of waveform X sample with period dt, using a window DT. 

Nwin = round(DT/dt);
if rem(Nwin,2)==0, Nwin=Nwin-1; end
[Nsam, Ncol] = size(X);
PS = [];
for icol=1:Ncol,
    RM = runmax(X(:,icol),Nwin);
    isteady = find(RM(1:Nsam-Nwin+1)==RM(Nwin:Nsam));
    PS = [PS; RM(isteady)];
end
hist(PS, Nbin);
xlabel('peak potential (mV)');
ylabel('# per bin');












