function ee = exist(dum,nam);
% Probetubecalib/exist - true of named probetubecalib object exists
%   exist(Probetubecalib(), 'Foo') returns True if Probetubecalib object
%   named Foo exists. Existsence is checked by checking whther its
%   datafiles are present on this computer.
%
%   NOTE 
%   The first input arg is a dummy for Methodizing this function.
%
%   See also Probetubecalib/fullfilename, Methodizing.

ee = logical(exist(fullfilename(dum,nam),'file'));





