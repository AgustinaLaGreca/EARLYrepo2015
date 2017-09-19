function B = JLbinint2(Uidx, doPlot, ids);
% JLbinint2 - analyze bin. interaction from I/C/B data. 
%   JLbinint(Uidx, doPlot)
%   JLbinint(Uidx, doPlot, [i_I i_C i_B]); % select ipsi/contra/binaural dataset

[doPlot, ids] = arginDefaults('doPlot/idx', nargout<1, []);

% identify partners
db = JLdbase(Uidx);
DB = JLdbase;
DB = DB([DB.icell_run]==db.icell_run); % select data from same cell 
DB = DB([DB.SPL]==db.SPL & [DB.freq]==db.freq); % same SPL & stim freq
db_i = DB([DB.chan]=='I');
db_c = DB([DB.chan]=='C');
db_b = DB([DB.chan]=='B');


if isempty(db_i),
    error('No Ipsi data available for triplet analysis.');
elseif isempty(db_c),
    error('No Contra data available for triplet analysis.');
elseif isempty(db_b),
    error('No Binaural data available for triplet analysis.');
end

if isempty(ids),
    if numel(db_i)>1,
        warning(['Multiple (' num2str(numel(db_i)) ') Ipsi data available; taking first one.']);
    elseif numel(db_c)>1,
        warning(['Multiple (' num2str(numel(db_c)) ') Contra data available; taking first one.']);
    elseif numel(db_b)>1,
        warning(['Multiple (' num2str(numel(db_b)) ') Bin data available; taking first one.']);
    end
    ids = [1 1 1];
end
[db_i, db_c, db_b] = deal(db_i(ids(1)), db_c(ids(2)), db_b(ids(3)));
Di = JLdatastruct(db_i.UniqueRecordingIndex);
Dc = JLdatastruct(db_c.UniqueRecordingIndex);
Db = JLdatastruct(db_b.UniqueRecordingIndex);

B = JLbinint(Di, Dc, Db, doPlot);















