function [D, Er] = locate(dum, ExpName, doRelax); 
% Experiment/locate - locate Experiment files.
%    D = locate(experiment, 'RG09211') returns the directory where the data
%    files of experiment RG09211 are located. An error is thrown when the
%    requested experiment is not found. Experiment files reside in the
%    <expName> subdir of the <Experimenter> subdir of Datadir.
%
%    locate(experiment, 'RG09211', 1) suppresses the experiment-not-found
%    error and returns '' instead.
%
%    See also Datadir, Experiment, DSread, Experiment/listall,
%    experiment/parenfolder. 

if nargin<3, doRelax=0; end
if ~isvarname(ExpName), ExpName = 'q____________________q'; end; % avoid unnecessary trouble with ExpName '' etc

Er = '';

D = folder(dummy(dum, ExpName));

if isempty(D) && ~doRelax,
    msg = ['Experiment ''' ExpName ''' not found in any subdir of ''' parentfolder(dum)  '''.'];
    title = 'File Error';
    errorh = errordlg(msg,title);
    uiwait(errorh);
    %error(['Experiment ''' ExpName ''' not found in any subdir of ''' parentfolder(dum)  '''.']);
end








