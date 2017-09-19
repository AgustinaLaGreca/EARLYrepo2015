BW = 1000; % Hz
CF = 600:1000:9600;
CF = 600;
dt = 10; % us sample period
StimDir = 'C:\MATLAB6p1\work\marcel\stim';
Dur = 5000; % ms dur
%----NL series-----


% Amp = 120;
% for ii=1:length(CF),
%     [w, S] = multisine(dt,Amp, Dur, CF(ii), BW, 5);
%     FN = fullfile(StimDir, ['marcel_NL' dec2base(ii,10,3)])
%     save([FN '.usr'], 'w', '-ascii');
%     save([FN '.param'], 'S', '-mat');
% end
    
    
%---DC series----
Amp = 20;
Dur = 10*Dur;
for ii=1:length(CF),
    [w, S] = multisine(dt,Amp, Dur, CF(ii), BW, 5);
    FN = fullfile(StimDir, ['marcel_DC' dec2base(ii,10,3)])
    save([FN '.usr'], 'w', '-ascii');
    save([FN '.param'], 'S', '-mat');
end






