function CH = ChunkClick(CH, figh);
% ChunkClick - list chunk that is clicked 
%   ChunkClick(CH) prompts the user to click on a waveform graph, then
%   displays the chunk(s) in array CH that span(s) the time poijnt clicked.
%   CH must have t0 and t1 fields.
%
%  ChunkClick(CH, figh) specifies the graphics handle of the figure.
%  Default is gcf.
%
%   See also DeriStats, ChunkPlot.

if nargin<2, 
    figh = gcf;
end

figure(figh);
[t, dum] = ginput(1);

if isempty(t), CH=CH(1:0); return; end
CH = CH(betwixt(t,[CH.t0],[CH.t1]));
if isempty(CH),
    display('No chunk here');
end

for ic=1:numel(CH),
    disp('=========');
    disp(CH);
    disp('=========');
end







