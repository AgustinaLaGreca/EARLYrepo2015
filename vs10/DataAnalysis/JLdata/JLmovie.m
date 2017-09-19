function B = JLmovie;
% JLmovie - create beat movie

global Jb214_6 Yi Yc Yb

if isempty(Jb214_6),
    Jb214_6 = local_read('RG10214', 6);
end

if isempty(Yi),
    [Yi, Yc, Yb] = JL_Poster(4.2);
end
[TTb, Ri, Rc, Rb] = local_interp(Yi,Yc,Yb);

%Mdir = 'C:\D_Drive\doc\Meetings and Talks\ARO\ARO-2011\beatmovie';
Mdir = '\\Ee1285a\MSO$';
Mfile = fullfile(Mdir, 'Beating4.avi');
MPfile = fullfile(Mdir, 'Beating4');

[Ri, Rc, Rb] = deal([Ri; Ri], [Rc; Rc], [Rb; Rb]);
TTb = timeaxis(Ri, diff(TTb(1:2)));

hfig=figure;
set(hfig,'units', 'normalized', 'position', [0.28 0.517 0.483 0.354])
set(hfig,'DoubleBuffer','on');
set(gca,'xlim',[-80 80],'ylim',[-80 80],...
    'NextPlot','replace','Visible','off');
mov = avifile(Mfile, 'COMPRESSION', 'None', 'quality', 100, 'FPS', 15);
%mov = avifile(Mfile, 'quality', 100, 'FPS', 15);
Tshift = linspace(0, 250, 400/4);
Twin = [0 10];
set(gca, 'color', [0.8 0.9 0.8]);
drawnow;
for ii=1:numel(Tshift),
    local_plot(TTb, Ri, Rc, Rb, Twin+Tshift(ii));
    F(ii) = getframe(gca);
    mov = addframe(mov,F(ii));
end
mov = close(mov);
Nrep=5; Psearch = 0; Bsearch = 0; RefFrame = 1;
%mpgwrite(F, get(hfig, 'Colormap'), MPfile, [Nrep, Psearch, Bsearch, RefFrame]);
mpgwrite(F, get(hfig, 'Colormap'), MPfile);


function local_plot(TTb, Ri, Rc, Rb, Twin);
dt = diff(TTb(1:2));
Nsam = round(diff(Twin)/dt);
isam = round(Twin(1)/dt)+(1:Nsam);
Tim = TTb(1:Nsam);
ylim([-0.5 2.5]);
plot(Tim, Ri(isam), 'color', [0 0.4 1], 'linewidth', 2.5);
set(gca, 'color', [0 0 0],'xtick', [], 'ytick', []);
ylim([-0.5 2.5]);
xplot(Tim, Rc(isam), 'color', [1 0 0.2], 'linewidth',2);
Yoffset = 1.15;
xplot(Tim, Rb(isam)+Yoffset, 'color', 0.7*[1 1 1], 'linewidth', 1.5);
xplot(Tim, Ri(isam)+Rc(isam)+Yoffset, 'color', [0 1 0.3], 'linewidth',2.5);
set(gca);
xlim([min(Tim) max(Tim)]);
set(gca, 'color', [0 0 0])
%==========================================
function Jb = local_read(ExpID, icell);
if isequal({'RG10214', 6}, {ExpID, icell}),
    JLreadBeats(ExpID, icell);   % 200 Hz & up, 100-Hz steps
    Jb = struct('B80', JbB_80dB, 'I80', JbI_80dB, 'C80', JbC_80dB, ...
        'B50', JbB_50dB, 'I50', JbI_50dB, 'C50', JbC_50dB);
elseif isequal({'RG10214', 4}, {ExpID, icell}),
    JLreadBeats(ExpID, icell);
    Jb = struct('B60', JbB_60dB, 'I60', JbI_60dB, 'C60', JbC_60dB);
else,
    error('don''t know what to read')
end

function [TTb, Ri, Rc, Rb] = local_interp(Yi,Yc,Yb);
Ri = repmat(Yi.ipsi_MeanWave-mean(Yi.ipsi_MeanWave), Yi.Freq1/4,1);
Rc = repmat(Yc.contra_MeanWave-mean(Yc.contra_MeanWave), Yc.Freq2/4,1);
Rb = Yb.bin_MeanWave-mean(Yb.bin_MeanWave);
TTi = Xaxis(Ri, Yi.ipsi_dt);
TTc = Xaxis(Rc, Yc.contra_dt);
TTb = Xaxis(Rb, Yb.bin_dt);
Ri = interp1(TTi, Ri, TTb);
Rc = interp1(TTc, Rc, TTb);












