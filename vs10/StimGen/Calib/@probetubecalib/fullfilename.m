function E = fullfilename(dum,nam);
% Probetubecalib/fullfilename - full file name of probe tube data file
%   fullfilename(Probetubecalib(), 'Foo') returns the full file name 
%   of probetibecalibe datafile Foo.
%
%   NOTE 
%   The first input arg is a dummy for Methodizing this function.
%
%   See also Probetubecalib/Datadir, Probetubecalib/fileExtension, 
%            Probetubecalib/exist, Methodizing.
E = fullfile(datadir(dum), [nam, fileExtension(dum)]);




