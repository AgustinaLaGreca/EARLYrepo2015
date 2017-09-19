function ConvertRueData(Nsam);
% ConvertRueData - read Rue's AM data from text files and save them as MAT

if isequal('Ee1285a', CompuName),
    Ddir = 'D:\data\RueData'
else,
    error('where are the data''s?');
end

qq = dir([Ddir '\*.txt']);
for ii=1:numel(qq),
    fn = qq(ii).name
    nam = strtok(fn,'.');
    nam = Words2cell(nam,'_');
    idata = 1+str2num(nam{3});
    nam = [nam{1} '_' nam{2}];
    wav = textread(fullfile(Ddir,fn),'%f', 'delimiter','\t');
    save([nam '_' num2str(idata) '.mat'], 'wav', '-mat');
end




    



