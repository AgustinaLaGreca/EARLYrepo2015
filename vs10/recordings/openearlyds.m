function DS = openearlyDS(FN);
%   OpenEarlyDS - look at Early data from file.
%      OpenEarlyDS('Foo.EarlyDS') displays the data n file Foo.EarlyDS.
%      
%   (experimental version)

[Path, Nam] = fileparts(FN);
[dum, Expname] = fileparts(Path);
iRec = str2num(Nam(end-4:end));
DS = DSread(Expname, iRec);
curds(DS);
dotraster(DS);



