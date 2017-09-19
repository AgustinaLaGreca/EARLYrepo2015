function [h, Lh, Bh]=plotpanes(Nwin, Arat, figh);
% plotpanes - divide plot window in N subplots or a few more
%    h=plotpanes(Nwin) finds integers N,M such that N/M is close to
%    3/4, and N*M>=Nwin, and creates subplots k=1:Nwin by invoking
%    h(k)=subplot(N,M,k).
%
%    plotpanes(Nwin, Arat) uses "aspect ratio" Arat instead of 3/4.
%
%    plotpanes(Nwin, Arat, figh) operates on the figure with handle figh.
%    Default is gcf.
%
%    plotpanes(Nwin, 0, figh) uses the height/width ratio of the figure as
%    Arat.
%
%    [h, Lh, Bh] = plotpanes(...) also returns the handles Lh of the N axes
%    at the left edge of the figure, and the handles Bh of the M axes that
%    have no lower neighbors. Typically, these axes are the ones to supply
%    with Y-labels and X-labels, respectively.
%    
%    See also SUBPLOT, Xlabels, Ylabels.

if nargin<2, Arat=3/4; end
if nargin<3, figh=gcf; end

if isequal(0,Arat),
    Fpos = get(figh,'position');
    Arat = Fpos(4)/Fpos(3);
end

N = ceil(sqrt(Nwin*Arat));
M = ceil(sqrt(Nwin/Arat));
if (N-1)*(M-1)>=Nwin,
    N = N-1; M = M-1;
elseif (N-1)*M>=Nwin && (~(N*(M-1)>=Nwin) || (N-1)*M<=N*(M-1) ) ,
    N = N-1;
elseif N*(M-1)>=Nwin,
    M = M-1;
end
        
%N,M,Nw=N*M
for ii=1:Nwin,
%     figure(figh);
    h(ii)= subplot(N,M,ii);
end
        
%N,M,Nwin, NM=N*M
Lh = h(1:M:end);
Bh = [h(end-M+(1:M))];
set(h(end), 'visible', 'off');








