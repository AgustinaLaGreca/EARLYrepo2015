function obj = struct2obj(S, Class);
% struct2obj - convert struct to object
%   struct2obj(S, 'myclass') converts struct S to object having requested
%   class. The fieldnames of S must match those of the object class.

obj = eval(Class); % void object
if ~isequal(fieldnames(obj), fieldnames(S)),
    error(['filenames of struct S do not match those of class ' Class]);
end

CirDir = cd;
cd(tempdir);
odir = ['@' Class];
if ~exist(odir, 'dir'),
    mkdir(odir);
end
FN = fullfile(tempdir, odir, 's2o.m');
if exist(FN,'file'),
    delete(FN);
end
textwrite(FN, 'function obj = s02(S, obj)');
textwrite(FN, 'FNS = fieldnames(S);');
textwrite(FN, 'for ii=1:numel(FNS),');
textwrite(FN, '   fn = FNS{ii};');
textwrite(FN, '   obj.(fn) = S.(fn);');
textwrite(FN, 'end');
rehash
obj = s2o(S, obj);
cd(CirDir);




