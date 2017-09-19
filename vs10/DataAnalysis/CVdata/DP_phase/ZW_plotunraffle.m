function ZW_plotunraffle(R, index, axhan)
%ZW_plotunraffle -make a plot of ZW_unraffleDPs.m's output
%
%ZW_plotunraffle(R, index) takes struct R from ZW_unraffleDPs and generates
%   a contourplot of the data in within the equal-Fdp band indicated by
%   index. The plot is created in a new figure.
%   Index may be an array, in which case each contourplot is generated in its
%   own new figure.
%
%ZW_plotunraffle(R,index,axhan) allows specification of an axes handle in
%   which to make the plot. Depending on the axes settings, the plot is
%   either added or replaces the current plot. Realise that when index is
%   an array, all the contours are plotted to this one axes.
%
% A contourplot is only generated if 6 or more significant DPs occur within
% the equal-Fdp band.
%
% See also ZW_unraffleDPs.

minDP = 6;
axhan = arginDefaults('axhan',[]);

iWrng = find(~ismember(index,1:R.Nband));
if ~isempty(iWrng),
    error(['Specified indices ' vector2str(index(iWrng)) ' are outside possible bands (1...' int2str(R.Nband) ')']);
end

%--- use recursion for multiple indices ---
if numel(index)>1,
    for ii = 1:numel(index),
        ZW_plotunraffle(R, index(ii),axhan);
    end
    return
end

%--- single index from here ---
if isempty(axhan), figure; axhan = axes; end

wght = R.wght{index}; %weight factors
iKp = wght>0;
%...if not enough DPs, skip plotting...
if sum(iKp)<minDP,
    if sum(iKp)==1,
        th(1) = text(1,1.2,'Only 1 DP is significant','parent',axhan);
    else,
        th(1) = text(1,1.2,['Only ' int2str(sum(iKp)) ' DPs are significant'],'parent',axhan);
    end
    th(2) = text(1,.8,'  No plot is generated','parent',axhan);
    set(th,'color','r','fontsize',18);
    axis(axhan, [0 4 0 2]);
    set(axhan,'visible','off');
    return
end

F = R.f{index}(iKp)/1e3; %f-freq; in kHz
G = R.g{index}(iKp)/1e3; %g-freq; in kHz
H = R.h{index}(iKp)/1e3; %h-freq; in kHz
data = R.data{index}(iKp); %the data

FG = F-G;
FG2H = F+G+2*H;
Fdp = F+G-H;

Ninterp = 20; 
[X, Y] = meshgrid(linspace(min(FG2H), max(FG2H),Ninterp), linspace(min(FG), max(FG),Ninterp));
% F=triscatteredInterp(FG2H,FG,data,'linear');
% Z = F(X,Y);
Z = griddata(FG2H,FG,data,X,Y);

cntrs = -1:.05:1;
clrmp = [circcolormap; circcolormap];
if isequal(R.datatype,'amplitude'),
    minZ = min(min(Z));
    maxZ = max(max(Z));
    cntrs = floor(minZ):3:ceil(maxZ);
    clrmp = colormap(axhan);
end

colormap(axhan,clrmp);
contourf(axhan, X,Y,Z,cntrs);
xplot(FG2H, FG,'ko','markersize',4,'markerfacecolor','w');
axis(axhan,[floor(min(FG2H)) ceil(max(FG2H)) floor(min(FG)) ceil(max(FG))]+[-1 1 -1 1]);
xlabel(axhan,'f+g+2h [kHz]');
ylabel(axhan,'f-g [kHz]');
title(axhan,[R.datatype '; equal-F_{dp} band #' int2str(index) '; (Fdp= ' num2str(round(mean(Fdp)*1e2)/1e2) ' kHz)']);
colorbar('peer',axhan);
axis(axhan,'square');
if ~isequal(R.datatype,'amplitude'), caxis([-1 1]); end










