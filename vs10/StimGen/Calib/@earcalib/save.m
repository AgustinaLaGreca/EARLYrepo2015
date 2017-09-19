function FFN = save(EC);
% Earcalib/save - save Earcalib object to file.
%    FFN = save(EC) saves T to file Foo.earcalib in the folder datadir(EC). 
%    The filename is equal to that of the experiment to which EC belongs,
%    together with a postfix like _001, _002, .. that allows multiple
%    calibrations to be measured within one experiment. The directory is
%    equal to datadir(EC). Save returns the full filename FFN.
%
%    In order to save EC, it must be filled with data (see Earcalib/isfilled) 
%    and belong to the current experiment - that is, iscurrent(EC.Experiment) 
%    must be true.
%
%    See also Earcalib/datadir, Earcalib/load.

if ~isfilled(EC),
    error('Earcalib object is not filled.');
elseif ~iscurrent(EC.Experiment),
    error('Cannot save Earcalib data whose experiment is not the current experiment.');
end

FileFilter = fullfile(datadir(EC), [name(EC.Experiment) '*' fileExtension(EC)]);
warning('off','MATLAB:warn_truncate_for_loop_index');
for ifile=1:inf,
    FFN = strrep(FileFilter, '*', ['_' dec2base(ifile,10,3)]);
    if ~exist(FFN, 'file'), 
        EC.iCalib = ifile; % index of calibration within experiment
        break;
    end
end
warning('on','MATLAB:warn_truncate_for_loop_index');
save(FFN,'EC','-mat');





