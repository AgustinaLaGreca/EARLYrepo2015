function AddMLsig;
% AddMLsig - add MLsig toolbox to EARLY path

%----high-level stuff
local_add2path('init');
% local_add2path('DA'); % D/A conversion
% local_add2path('Joystick'); % pointer stuff
% local_add2path('TDT'); % TDT drivers
% local_add2path('TDT\SYS3'); % friggin' sys 3 stuff
% local_add2path('2AFC'); % forced-choice procedures
% local_add2path('DataHandling'); % book keeping for experiments
% local_add2path('LesGrant'); % fringe experiments from Les' grant
%----utility functions
% local_add2path('util\tools');
% local_add2path('util\strfun');
% local_add2path('util\misc');
local_add2path('util\conversion');
local_add2path('util\errorfun');
%----lower-level stuff
local_add2path('MLsigObject'); % MLsig object definition
local_add2path('MLsigObject\tools');
local_add2path('MLsigObject\MLelfun'); % MLsig elem functions
local_add2path('MLsigObject\MLdatafun');% MLsig data functions
local_add2path('MLsigObject\MLsignal');% MLsig signal functions
local_add2path('MLsigObject\create');% MLsig signal generators
local_add2path('MLsigObject\demos');% audio demos

%---------locals-------------
function local_add2path(D);
addpath(fullfile(versiondir, 'MLsig', D), '-end');



