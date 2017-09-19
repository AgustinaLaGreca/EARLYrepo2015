function tshift = templot(Y, Threshold, doPlot);
% templot - aligned plot of action potential templates
%    Tshift=templot(Y, Threshold, doPlot), where Y is the return arg of 
%    xspikes, plots the spike templates while aligning their first 
%    crossing of threshold value Threshold. Return arg Tshift is the
%    category-wise time shift needed to align the templates.
%
%    If doPlot is false, plotting is supressed, but Tshift is returned.
%
%   See also xspikes.

if nargin<3, doPlot=1; end

Ntem = numel(Y.EventCount);
dt = diff(Y.t0snip(1:2));
for item=1:Ntem,
    tem = Y.Template{item}+Y.Ymean;
    NsamTem = numel(tem);
    imid = round((NsamTem+1)/2);
    ithr = find(tem>Threshold,1);
    if isempty(ithr), continue; end
    tshift(item) = 0.3 + dt*(ithr-imid); % time shift
    if doPlot,
        %xplot(Y.t0snip, tem, lico(item), 'linewidth', 2);
        xplot(Y.t0snip-tshift(item), tem, lico(item), 'linewidth', 2);
    end
end
xplot(xlim, [1 1]*Threshold, 'k:');









