function [Y, dt, Dur] = ReadRueData(ID, irep, Ncond, TimeWin);
% ReadRueData - read Rue's data from IBW file
%  [Y, dt, Dur] = ReadRueData(ID, irep, Ncond, TimeWin);
%    default Ncond = [] meaning all stim conditions.
%    TimeWin is time window in ms. Default = 450 ms = whole resp
% 
%  Examples
%       [Y, dt, Dur] = ReadRueData('p90605Bi',0, [], 450)
%    rep #0 read from p90605Bi_0_.ibw
%     Dur  =  104400 ms = 8 x 29 x 450 ms. (8 SPLs, 29 freqs)
%
%        [Y, dt, Dur] = ReadRueData('p90605B*',0, [], 450)
%    or, equivalently,
%        [Y, dt, Dur] = ReadRueData('p90605Bq',0, [], 450)
%    concatenates ipsi and contra responses
%    rep# read from p90605Bi_0_.ibw and p90605Bc_0_.ibw.



%  [Y, dt, Dur] = ReadRueData(ID,irep, Ncond, TimeWin);
%    default 

if nargin<3, Ncond=[]; end % default : all conditions
if nargin<4, TimeWin=450; end % ms analysis window per stimulation

if isempty(Ncond), Ncond=29*8; end

NsamPres = round(25*450); % # samples in each pres
Ntw = round(25*TimeWin); % # selected samples per pres

if Ntw>NsamPres, error('TimeWin may not exceed stim dur.'); end

if isequal('Ee1285a', CompuName) || isequal('SIUT', CompuName),
    Ddir = 'D:\data\RueData';
else,
    error('where are the data''s?');
end

Ext = '.mat'; % default extension
switch ID,
    case 'A',
        FN = '090305_A';
    case 'B',
        FN = '090305_B';
    case 'Z',
        FN = '090206_A';
    otherwise,
        FN = ID; % literal, e.g. p90605Bi
        ismat = 0;
        if ~isempty(strfind(ID,'FM')),
            Ddir = 'D:\Data\RueData\IBW\FM';
        else,
            Ddir = 'D:\Data\RueData\IBW\Audiograms';
        end
        Ext = '.ibw';
end

if ismat,
    FN = fullfile(Ddir, [FN '_' num2str(irep) Ext]);
    Y = load(FN);
    Y = Y.wav;
else, % IBW file; full file name given
    FN = fullfile(Ddir, [FN '_' num2str(irep) '_' Ext ]);
    if ismember('*', FN) || ismember('q', FN), % concatenate ipsi and contra
        FNi = strrep(strrep(FN,'*', 'i'), 'q', 'i');
        FNc = strrep(strrep(FN,'*', 'c'), 'q', 'c');
        Y(1) = IBWread(FNi);
        Y(2) = IBWread(FNc);
        Y = [Y(1).y; Y(2).y];
        Ncond = 2*Ncond;
    else,
        Y = IBWread(FN);
        Y = Y.y;
    end
end

Y = Y(1:(Ncond*NsamPres)); % restrict to requested # conditions
%NsamPres, Ncond ,NsamPres*Ncond , Nsam
Y = reshape(Y,NsamPres,Ncond);  %each col is one pres
Y = Y(1:Ntw,:); % select restricted time window
Y = Y(:);

dt = 1/25; % ms sample period
Dur = numel(Y)*dt;



    



