function ZWOAEstruct
%ZWOAESTRUCT -->>help ZWOAEstruct: explains different fields within data-structure S for ZWOAE
%
%The structure S which holds all info for a single recording of ZWOAE data
%holds the following fields:
%
% .EXPname        : name of the experiment, e.g. 'RG07046'
% .CALname        : name of calibration data file, e.g. 'speaker1.cal'
% .index          : integer identifying the record number within EXPname
% .date           : date of the original recording
% .time           : time of the original recording
%
% .fs             : sample frequency [kHz]
% .numaverage     : number of blocks that were averaged
% .micsensitivity : struct with fields that speccify the microphone
%                   sensitivity. Fields are:
%                   .dBSPL:      : A tone of "dBSPL" dB SPL gives a mic. output of "dBV" dBV
%                   .dBV         : output of microphone for a tone of "dBSPL" dB SPL
%
% .signal         : averaged signal, first periodic block is omitted to
%                   exclude onset phenomena
% .micgain        : gain setting for microphone signal
% .units          : structure with units for different saved parameters.
%                 Fields are:
%                    .fsunit      : unit for the sample frequency,         default = 'kHz'
%                    .frequnit    : unit for the frequencies (f1,f2, fdp), default = 'Hz'
%                    .Levelunit   : unit for the levels (L1,L2),           default = 'dB SPL'
%                    .phaseunit   : unit for the phases,                   default = 'rad'
%                    .signalunit  : unit for the averaged signal,          default = 'V'   
%                    .gainunit    : unit for microphone gain,              default = 'dB'
%
% .stimpars       : structure with "command line" parameters. Here, the
%                  amplitude of several freq. components may be very low (i.e. not played)
%                  Within .stimpars are the following fields:
%                    .N1          : array with on/off bit mask for all
%                                   F1's (1 is on, 0 is off)
%                    .N2          : array with on/off bit mask for all
%                                   F2's (1 is on, 0 is off)
%                    .F1          : array with all F1-frequencies
%                    .F2          : array with all F2-frequencies
%                    .L1          : array with all L1-levels (as parameter)
%                    .L2          : array with all L2-levels (as parameter)
%                    .phi1        : array with all starting phases for F1's
%                    .phi2        : array with all starting phases for F2's
%
% .setupinfo      : structure holding info about the recording setup,
%                   experimenter, and the animal. Fields are:
%                    .ear         : which ear? left ('l') or right ('r') ear               
%                    .numchan     : which D/A channels used? e.g. [1],[2],[1 2],
%                    .mictype     : what microphone? e.g. 'ER-10A'
%                    .drivertype  : what drivers/speakers? e.g. 'ER-1'
%                    .DAtype      : what D/A-A/D converter? e.g. 'RP2','RX6',
%                    .attentype   : if attenuators and for which channels,
%                                   e.g. [0], [0 0], [1 0], [1 1]
%                    .experimenter: who did the recordings? e.g. 'Bas'
%