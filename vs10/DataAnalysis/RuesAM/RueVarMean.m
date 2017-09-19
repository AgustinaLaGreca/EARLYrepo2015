function S = RueVarMean(Id);
%  RueVarMean(Id);

Nfile.A = 20; Nfile.B = 10; Nfile.Z=9;


Nseg = 10;
NsamSeg = 2610000/Nseg;
M = []; V = [];
for iseg=1:Nseg,
    iseg
    offset = (iseg-1)*NsamSeg;
    Seg = [];
    for ii=1:Nfile.(Id),
        s = ReadRueData(Id,ii);
        Seg = [Seg, s(offset+(1:NsamSeg))];
    end
    dsize(Seg);
    M = [M; mean(Seg,2)];
    V = [V; var(Seg,[],2)];
end

S = CollectInStruct(M,V);















;