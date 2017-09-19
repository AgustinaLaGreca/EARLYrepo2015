% Sequence playback on Sys 3 devices
%
% Sequenced play
%   seqplayinit           - initialize TDT device for sequenced play
%   seqplayupload         - Upload samples for sequenced play over TDT device
%   seqplaylist           - specify play list for sequenced play
%   seqplaygo             - trigger sequenced play over sys 3 device
%   seqplayhalt           - immediately halt sequenced play D/A
%   seqplaystatus         - current status of sequenced play over Sys3 Device
%   SeqplayDAwait         - wait for seqplay DA conversion to finish
%   SeqplayGrab           - set or get grabbing status of seqplay circuit
%
% Diagnostics & debugging
%   TimingCalibrateRX6    - calibrate timing between analog output & digital input of RX6
%   sing                  - audio test using Seqplay functionality
%   seqplayinfo           - current status of Seqplay variables
%   SeqPlayPlotStatus     - plot sequence play history (debug function)
%   seqplayplot           - plot buffer contents of seq play (debug function)
