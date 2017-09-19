function hax = deriplot(RW, AP, EPSP, Chunk),
% deriplot - plot deristat output
%    deriplot(RW, AP, EPSP, CH), where [CH, RW, AP, EPSP] are the output
%    args of deristat. Note the change of order! CH is optional. Deriplot
%    returns the handle to the axes.
%
%    See also deristat, deriplot.

if nargin<4, Chunk=[]; end

[sam, T] = RW.getsam();
dt = RW.dt;
D = readTKABF(RW); 
StimProps = XgetStimprops(D);
BaseLines = RW.MedSam+RW.IqrSam*[-1 1]*RW.IqrFactor(1)/2;
t2isam = @(t)(1+round(t/dt));
t2sam = @(t)sam(1+round(t/dt));

% plot waveform, highlighting Chunks, if any
figh = BasicChunkPlot(Chunk, dt, sam, BaseLines);
hax = findobj(figh,'type', 'axes');
% indicate stimulus bursts
Onset = StimProps.onset; Offset = StimProps.offset; Nstim = numel(Onset);
dashplot((Offset+Onset)/2, RW.MedSam+0*Onset, (Offset-Onset), 'color', [1 0.8 0.8], 'linewidth', 3);

if ~isempty(Chunk),
    tAPdown = [Chunk([Chunk.isAPdownFlank]).tmaxRate];
    xplot(tAPdown, t2sam(tAPdown), 'ko');
    chunkPlot(Chunk([Chunk.isAPdownFlank]), T, sam, 'k', 'linewidth', 2);
    %ChunkPlot(Chunk([Chunk.isEPSP]), T, sam, 'g');
    chunkPlot(Chunk([Chunk.isAPpreUpFlank]), T, sam, 'g', 'linewidth', 2);
    chunkPlot(Chunk([Chunk.isAPpostUpFlank]), T, sam, 'c', 'linewidth', 2);
end

APdown = struct('istart', num2cell(t2isam([AP.tpeak])), 'iend', num2cell(t2isam([AP.tdip])));
tAPdown = [AP.tmaxDownSlope];
xplot(tAPdown, t2sam(tAPdown), 'ko');
chunkPlot(APdown, T, sam, 'k', 'linewidth', 2);
APpostUp = struct('istart', num2cell(t2isam([AP.tdip])), 'iend', {AP.iend});
chunkPlot(APpostUp, T, sam, 'c', 'linewidth', 2);
APpreUp = struct('istart', num2cell(t2isam([AP.tInflex])), 'iend', num2cell(t2isam([AP.tpeak])));
APpreUp = APpreUp(~isnan([APpreUp.istart]) & ~isnan([APpreUp.iend]));
chunkPlot(APpreUp, T, sam, 'g', 'linewidth', 2);

chunkPlot(EPSP(~[EPSP.isAPpostUpFlank]), T, sam, '-', 'color', [0.1 0.6 0.1], 'linewidth', 2);
chunkPlot(EPSP([EPSP.isAPpostUpFlank]), T, sam, ':', 'color', [0.1 0.6 0.1], 'linewidth', 2);
xplot([EPSP.tpeak], [EPSP.Vpeak], '*', 'color', [0.1 0.6 0.1]);

% Cr = Chunk([Chunk.isrising]& (~[Chunk.isAPupFlank]) & ([Chunk.maxRate]<ThrEPSPrate));
% CrN = Cr([Cr.y0]<MedSam-2*IqrSam);
% ChunkPlot(CrN, T, sam-0.02*diff(ylim), 'k-');

fenceplot([AP.tInflex], ylim, 'color', 0.75*[1 1 1])
%xplot(tinflx, t2sam(tinflx), 'k+', 'markersize', 8);

Tstr = [D.ExpID ' <' D.RecID '> icond=' num2str(D.icond) ' --- ', ...
    num2sstr(StimProps.Fcar) ' Hz  ' num2str(StimProps.SPL) ' dB' ];
title(Tstr);
setGUIdata(gcf, 'StimProps', rmfield(StimProps, 'DS'));
TracePlotInterface(gcf); % navigation keys








