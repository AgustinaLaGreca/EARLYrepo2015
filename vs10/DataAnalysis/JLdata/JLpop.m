%=======JLpop=======
% population analysis of JLbeat data

% select 400-Hz, binaural data
DB400Hz_B = DB([DB.freq]==400 & [DB.chan]=='B');
struct2char(DB400Hz_B) % display list to screen
Jds = JLdatastruct(DB400Hz_B);
% select those w monaural partners
ifull = find([Jds.hasMonPartners]);
DB400Hz_B = DB400Hz_B(ifull);
JdsFull = Jds(ifull);

% save and load
save D:\processed_data\JL\JLpop DB400Hz_B JdsFull
load D:\processed_data\JL\JLpop DB400Hz_B JdsFull

FNS = fieldnames(db)';
FNS = [FNS(1:end-1) {'I_partners' 'C_partners' 'B_partners'}];
JdsX= structpart(JdsFull,FNS)
sp2sc = @(s)strrep(trimspace(num2str(s)),' ','_');
IP = cellfun(sp2sc, {JdsX.I_partners}, 'unif',0);
CP = cellfun(sp2sc, {JdsX.C_partners}, 'unif',0);
BP = cellfun(sp2sc, {JdsX.B_partners}, 'unif',0);
[JdsX.sI_partners] = deal(IP{:});
[JdsX.sC_partners] = deal(CP{:});
[JdsX.sB_partners] = deal(BP{:});

struct2char(rmfield(JdsX,{'I_partners' 'C_partners' 'B_partners' }))

%=======================================
%=======================================
Pdir = fullfile(processed_datadir, 'JL', 'JLpop');
[D3, P] = JLmonpartners(JLdbase); % complete triangles
D3P = structJoin(D3, P);
save(fullfile(Pdir, 'P_D3'), 'P', 'D3', 'D3P');
load(fullfile(Pdir, 'P_D3'));
D3_400Hz = D3([D3.freq]==400);
P_400Hz = P([D3.freq]==400);
D3P_400Hz = P([D3P.freq]==400);
% complete JLbinint; results are cached
%>>>> rmcache binint*    clears cache
for ii=1:numel(P); 
    ii, 
    JLbinint(P(ii).ipsi_partner(1), P(ii).contra_partner(1), D3(ii).UniqueRecordingIndex,0); 
end
% -----------------------------------

% -------JLbinint metrics
clear JbM 
for ii=1:numel(P),
    ii
    JbM(ii) = JLbinintMetrix(P(ii).ipsi_partner(1), P(ii).contra_partner(1), D3(ii).UniqueRecordingIndex); 
end
save(fullfile(Pdir, 'JbM'), 'JbM');
load(fullfile(Pdir, 'JbM'), 'JbM');
% selection for following analyses
jm = JbM;
jm = JbM([JbM.freq]==400);
% Pk-Pk of monaural recs vs Pk-Pk of mon contrib to binaural recs
Si = JbMscatter(jm, 'ipsi_P2P_ipsi', 'bin_P2P_ipsi', 'b.');
Sc = JbMscatter(jm, 'contra_P2P_contra', 'bin_P2P_contra', 'r.');
xplot(xlim, xlim, 'k');
legend({['R_{ipsi}=' num2str(Si.corr)], ['R_{contra}=' num2str(Sc.corr)]}, 'location', 'best');
xlabel('Monaural Pk-Pk'); ylabel('Binaural Pk-Pk');
% mon vs bin shapes
Si = JbMscatter(jm, 'Rmax_ipsi', 'Rmax_contra', 'k.');
xplot(xlim, xlim, 'k');
xlabel('Ipsi corr(mon,bin)'); ylabel('Ipsi corr(mon,bin)');


% #peaks vs Pk-Pk
plot([JbM.ipsi_P2P_ipsi], [JbM.ipsi_Npeak_ipsi], '.')
xplot([JbM.contra_P2P_contra], [JbM.contra_Npeak_contra], 'r.')


% view all bin interaction analyses of 400-Hz triples
d3p = D3P([D3P.freq]==400); % restrict triple database
for ii=1:numel(d3p),
    JLbinint(d3p(ii).ipsi_partner(1), d3p(ii).contra_partner(1), d3p(ii).UniqueRecordingIndex);
    pause;
    aa;
end

% view all bin interaction analyses of one cell
d3p = D3P([D3P.icell_run]==33); % restrict triple database
for ii=1:numel(d3p),
    JLbinint(d3p(ii).ipsi_partner(1), d3p(ii).contra_partner(1), d3p(ii).UniqueRecordingIndex);
    pause;
    aa;
end


