function CC = RueCorr(Id1, Id2, Ncond, TimeWin);
%  RueCorr(Id1, Id2);

if nargin<3, Ncond=29*8; end
if nargin<4, TimeWin=450; end % ms analysis window per stimulation

Nfile.A = 20; Nfile.B = 10; Nfile.Z=9;

N1 = Nfile.(Id1); N2 = Nfile.(Id2);
sameSet = isequal(Id1,Id2);

more off;
CC = [];
for i1=1:N1,
    disp([num2str(i1) '/' num2str(N1)]);
    W1 = ReadRueData(Id1,i1, Ncond, TimeWin); 
    if sameSet, i2start=i1+1; else, i2start=1; end;
    for i2=i2start:N2,
        disp(['   ' num2str(i2) '/' num2str(N2)]);
        W2 = ReadRueData(Id2,i2, Ncond, TimeWin);
        cc = corrcoef(W1,W2);
        CC = [CC, cc(1,2)];
    end
end













