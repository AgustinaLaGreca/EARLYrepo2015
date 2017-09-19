function SPMcycleHist(MM, Nbin, PHshift, XLpos, Xsh);
% function SPMcycleHist - cycle histogram from spike metrics output
%   SPMcycleHist(M);
if nargin<3, PHshift=[]; end
if nargin<4, XLpos = 0.1; end % norm X pos to put text
if nargin<5, Xsh = 0; end

[PHe, PHa] = deal([]);
for im = 1:numel(MM),
    M = MM(im);
    T = 1e3/M.Fcar;
    du= M.duringStim;
    T_ep = M.itvPrevOnset(du);
    T_ap = denan(T_ep + M.itvEPSP_AP(du));
    [dum, Phep(im)] = vectorstrength(1e-3*T_ep, M.Fcar);
    [dum, Phap(im)] = vectorstrength(1e-3*T_ap, M.Fcar);
    if isempty(PHshift),
        dPHe = -Phep(im)+1.4;
        dPHa = -Phap(im)+1.4;
    else,
        [dPHe, dPHa] = deal(PHshift);
    end
    ph_ep = rem((T_ep/T) + dPHe, 1);
    ph_ap = rem((T_ap/T) + dPHa, 1);
    PHe = [PHe; ph_ep(:)];
    PHa = [PHa; ph_ap(:)];
end

Rep = vectorstrength(PHe, 1);
Rap = vectorstrength(PHa, 1);

PHe = rem(PHe + Xsh(1), 1);
PHa = rem(PHa + Xsh(end), 1);


figure
dphi = 1/Nbin;
set(gcf,'units', 'normalized', 'position', [0.476 0.302 0.242 0.569]);
subplot(2,1,1);
hist(PHe, linspace(dphi/2, 1-dphi/2,Nbin));
xlim([0 1]);
text(XLpos(1), 0.9*max(ylim), ['All events: R=' num2str(Rep,2)]);
ylabel('# events');

subplot(2,1,2);
hist(PHa, linspace(dphi/2, 1-dphi/2,Nbin));
xlim([0 1]);
text(XLpos(end), 0.9*max(ylim), ['APs only: R=' num2str(Rap,2)]);
ylabel('# events');
xlabel('Phase (cycles)');




