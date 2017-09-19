function S = JLcellProps(D);
% JLcellProps - determine & retrieve cell properties
%    JLcellProps(JLdbase);
%    CP = JLcellProps;


Props = {'ID' 'MHcomment' '-Stim' 'Nrec' 'Nseries' 'SPLi' 'SPLc' 'SPLb' 'SPLicb' '-Rates' 'Rates'};
if nargin==1, % compile
    Icell_run = unique([D.icell_run]);
    for icell=1:numel(Icell_run),
        icell
        ihit = find([D.icell_run]==Icell_run(icell));
        d = D(ihit(1));
        for iprop=1:numel(Props),
            prop = Props{iprop};
            if isvarname(prop),
                S(icell).(prop) = feval(['local_' prop], D(ihit));
            else, % header
                S(icell).([prop(2:end) '__________']) = '_________________';
            end
        end
    end
end
S = tree2hedge(S); % "flatten" the struct, i.e., unpack its struct-valued fields
% ============================================

function Y = local_ID(D);
Y = structpart(D(1),{'iexp' 'icell' 'icell_run'});

function Y = local_MHcomment(D);
[Y, CFN, CP] = getcache([mfilename '_local_MHcomment'], [D.UniqueRecordingIndex]);
if ~isempty(Y), return; end; clear Y;
S = JLdatastruct(D(1));
Y = S.MHcomment;
putcache(CFN, 200, CP, Y);

function Y = local_Nrec(D);
Y = numel(D);

function Y = local_Nseries(D);
Y = numel(unique([D.iseries_run]));

function Y = local_SPLi(D);
D = D([D.chan]=='I');
Y = unique([D.SPL]);

function Y = local_SPLc(D);
D = D([D.chan]=='C');
Y = unique([D.SPL]);

function Y = local_SPLb(D);
D = D([D.chan]=='B');
Y = unique([D.SPL]);

function Y = local_SPLicb(D);
Di = D([D.chan]=='I');
Si = unique([Di.SPL]);
Dc = D([D.chan]=='C');
Sc = unique([Dc.SPL]);
Db = D([D.chan]=='B');
Sb = unique([Db.SPL]);
Y = intersect(Si,Sc);
Y = intersect(Y, Sb);

function Y = local_Rates(DB);
[Y, CFN, CP] = getcache([mfilename '_local_Rates'], [DB.UniqueRecordingIndex]);
if ~isempty(Y), return; end; clear Y;
[Nspike_spont, Nspike_driven, Dur_spont, Dur_driven] = deal(0); % no spikes, no duration yet
for ii=1:numel(DB),
    S = JLdatastruct(DB(ii)); % standardized recording info
    [D, DS, L]=readTKABF(S);
    dt = D.dt_ms; % sample period in ms
    rec = D.AD(1).samples; % entire recording, including pre- and post-stim parts
    Time = Xaxis(rec, dt);
    cutWin = [S.APwindow_start S.APwindow_end];
    [dum, SPTraw] = APtruncate2(rec, [S.APthrSlope, S.CutoutThrSlope], dt, cutWin, 1); 
    Nspike_spont = Nspike_spont + sum((SPTraw<S.preDur) | (SPTraw>S.preDur+S.burstDur+20));
    Nspike_driven = Nspike_driven + sum(betwixt(SPTraw, S.preDur, S.preDur+S.burstDur));
    Dur_spont = Dur_spont + S.preDur + max(Time)- (S.preDur+S.burstDur+20); % add spont duration of this rec
    Dur_driven = Dur_driven + S.burstDur; % add spont duration of this rec
end
Spont_rate = 1e3*Nspike_spont/Dur_spont; % spont spike rate (sp/s)
[Dur_spont, Dur_driven] = DealElements(round([Dur_spont, Dur_driven]));
Y = CollectInStruct(Nspike_spont, Nspike_driven, Dur_spont, Dur_driven, Spont_rate);
putcache(CFN, 200, CP, Y);


