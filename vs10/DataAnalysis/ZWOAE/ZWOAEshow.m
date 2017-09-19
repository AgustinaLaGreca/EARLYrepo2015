function D = ZWOAEshow(expID, recID, DPtype, ah)
% ZWOAEshow -plot amplitude spectrum of ZWOAE with approp. markers
%    D = ZWOAEshow(expID, recID, DPtype, ah) plots the amplitude spectra for the
%    specified recordings.
%
%    INPUT variables:
%
%    expID: numeric ID for the experiment/gerbil
%    recID: numeric ID for the recordings (may be an array)
%    DPtype: string specifying which DPs to consider, default = 'all'
%            options: 'near', 'far','supall','all'
%    ah:    handle to axes in which to generate the plot, If omitted, new
%           figure and axes are created.
%
%   OUTPUT variable:
%
%   D = struct as returned by ZWOAEimport. D is an array of structs when
%   recID specifies multiple recordings.
%
%   See also ZWOAEimport, ZWOAEplotstruct, ZWOAEdominance, SpecComp

defDPtype = 'all';
defah = '';

if nargin < 4, ah = defah; end
if nargin < 3 | isempty(DPtype), DPtype = defDPtype; end

%handle multiple datasets recursively
if numel(recID) > 1,
    for ii = 1:numel(recID),
        D(ii) = ZWOAEshow(expID, recID(ii));
    end
    return;
else
    D = ZWOAEimport(expID, recID);
end

%======== handle single data set from here =============
PS = ZWOAEplotStruct; %struct of structs with plot info

Nzw = D.Nzwuis;
Fprim = [D.Fzwuis(:); D.Fsingle];
% get correct DP matrix, compute DP freqs & get spectrum @ these freqs
[SM, SF] = ZWOAEmatrices(Nzw);
M = SM.(['M' DPtype]); % select requested DP type
Fdp = M*Fprim; %calculate DP frequencies
Fdp = abs(Fdp); %"fold" negative freqs into positive part of spectrum

FreqAx = linspace(0, 3*max(Fprim), 1000);
%Rspec = polyval(D.Pnf, FreqAx); %noise floor
Rspec = ZWOAEgetNoise(D, FreqAx);
Sprim = SpecComp(D.df, D.MG, Fprim); %amplitudes for primaries
Sdp = SpecComp(D.df, D.MG, Fdp); %amplitudes for DPs

%Look at how dominant ZWOAEs are
[SNRa, SNRb] = ZWOAEdominance(D.MG, D.df, D.Fzwuis,D.Fsingle, DPtype);

%make string for title to know what was plotted
FN = ZWOAEfilename(D.ExpID, D.recID);
Tstr = [int2str(D.ExpID) ' (# ' int2str(D.recID) ') ; '];
Tstr = [Tstr 'RecType: '  D.RecType(1:2)];
Tstr = [Tstr '; Z = ' num2str(round(D.Fz_mean*10)/10)];
Tstr = [Tstr ' (' int2str(D.Nzwuis) '); S = ' num2str(round(D.Fsingle*10)/10)];
Tstr = [Tstr '; Lzw = ' vector2str(D.Lzwuis)];
Tstr = [Tstr '; Ls = ' vector2str(D.Lsingle)];
Tstr = [Tstr '; DPtype = ' DPtype];

%create axes when necessary
if isempty(ah),
    figure;
else
    axes(ah);
end

%plot...
xdplot(D.df, D.MG, PS.Spec); %spectrum
xplot(FreqAx, Rspec, PS.Floor); %noise floor
xplot(Fprim, Sprim, PS.Prim); %mark primaries
xplot(Fdp, Sdp, PS.(['DP' DPtype])); %mark DPs

%determine necessary range on x-axis
tmp = SM.('Mall'); Fall = tmp*Fprim;
[xl(1) xl(2)] = minmax(Fall); xl = xl.*[.9 1.1];
xlim(xl);

%do similar thing for y-axis
yl(1) = min(Rspec)-20; yl(2) = max(Sprim+10);
ylim(yl);

%place dominance numbers within axes
text(xl(1)+0.65*diff(xl), yl(1)+0.9*diff(yl),['SNR all: ' num2str(SNRa),' dB']);
text(xl(1)+0.65*diff(xl), yl(1)+0.83*diff(yl),['SNR limit: ' num2str(SNRb),' dB']);

%add labels on axes
yh = ylabel('Level [dB SPL]');
xlabel('Freq [kHz]');

yp = get(yh,'position'); %get position of ylabel
tmp = get(gca,'title'); pp = get(tmp,'position'); %get position of title

text(yp(1)*.95, pp(2)*1.05, pp(3), Tstr,'Tag','TitleText'); %use it to place text
set(gca,'Layer','top'); box on;

%========obsolete attempt to mark all kinds of DPs==========
% return
% % plot markers on all 2nd & 3rd order DP's
% Fprim = [D.Fzwuis(:); D.Fsingle; D.Fsup];
% Fprim = Fprim(Fprim>0);
% qq=moam(numel(Fprim), 3);
% 
% [SM, SF] = ZWOAEmatrices(Nzw);
% M = SM.(['M' DPtype]); % select requested DP type
% 
% if ~isequal(D.Fsup, 0),
%     M = [M, zeros(size(M,1),1)];
% end
% qq = setdiff(qq, M,'rows');
% 
% tmp = sum(abs(qq),2);
% 
% Fall = qq*Fprim;
% Sall = SpecComp(D.df, D.MG, Fall);
% 
% tmp(Fall<=0)=[];
% Sall(Fall<=0)=[];
% Fall(Fall<=0)=[];
% 
% i2nd = find(tmp==2); %indices for 2nd-order
% i3rd = find(tmp==3); %indices for 3rd-order
% i4th = find(tmp==4); %indices for 4th-order
% i5th = find(tmp==5); %indices for 5th-order
% 
% xplot(abs(Fall(i2nd)), Sall(i2nd),'s','color',.8*[1 0 1],'markersize',5);
% xplot(abs(Fall(i3rd)), Sall(i3rd),'v','color',.7*[0 1 0],'markersize',5);




