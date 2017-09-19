function ChunkScatter(C, fn1, fn2, plotArgs, ChunkName, hAx);
% ChunkScatter - right-clickable scatter plot of chunk array
%     ChunkScatter(C, fn1, fn2, plotArgs, ChunkName, hAx) plots [C.(fn1)] 
%     against [C.(fnw)] using plotArgs as trailing arg to PLOT.
%     The markers of the plot are right-clickable. Selecting the "view"
%     label of the context menu will show 30 ms centered around the chunk
%     in the plot in axes system with handle hAx. Defafault hAx is gca.
%
%     C must be chunk struct array, containing at least fields t0, ..
%
%     See also deristats.

X = [C.(fn1)];
Y = [C.(fn2)];
hl = plot(X,Y,plotArgs);
IDpoints(hl, C, [C.t0], @(C,t)[ChunkName ' @ ' num2str(t) ' ms'], 'view', {@local_centerview});

%====================
function local_centerview(C, )



