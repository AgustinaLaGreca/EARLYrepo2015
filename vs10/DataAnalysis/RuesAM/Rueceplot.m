function Rueceplot(CC, minCorr, XX, YY, Cbar)
% Rueceplot - plot Rueceptive plot
%   Rueceplot(CC) plots Rue's receptive field from audiogram data.
%
%   Rueceplot(CC, minCorr) uses a range of minCorr to +1 for the mapping
%   the correlation values in CC to the color map. Default, minCorr = -1.
%
%   Rueceplot(CC, minCorr, XX, YY, Cbar) has additional logical inputs:
%      XX: wether to draw xlabel & xticklabels (Freq axis).
%      YY: wether to draw ylabel & yticklabels (SPL axis).
%    Cbar: whether to draw a colorbar.
%   By default, XX, YY, and Cbar are all true.
%
%   See also RueceptiveField, RueCross.

if nargin<2, minCorr=-1; end
if nargin<3, XX=1; end
if nargin<4, YY=1; end
if nargin<5, Cbar=1; end

CC = flipud(CC);
[NSPL, Nfreq] = size(CC);
CC = [[CC, ones(NSPL,1)]; minCorr*ones(1,Nfreq+1)]; % invisible outer elements of matrix do determine colormap range
pcolor(CC);
if XX,
    set(gca,'xtick', 1.5:5:29.5, 'xticklabel', {'1' '2' '4' '8' '16' '32'});
    xlabel('Frequency (kHz)');
else,
    set(gca,'xticklabel',{});
end
if YY,
    set(gca,'ytick', 1.5:9.5, 'yticklabel', {'0' '10' '20' '30' '40' '50' '60' '70' '80' });
    ylabel('Sound level (dB SPL)');
else,
    set(gca,'yticklabel',{});
end
if Cbar,
    colorbar;
end
figure(gcf);
colormap(RueCmap);




