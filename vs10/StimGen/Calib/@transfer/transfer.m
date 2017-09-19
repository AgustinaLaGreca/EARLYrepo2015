function T = transfer(Q_stim, Q_resp, Ref_stim, Ref_resp, dBref_stim, dBref_resp);
% transfer - constructor for Transfer objects
%   transfer(Q_stim, Q_resp, Ref_stim, Ref_resp, dBref_stim, dBref_resp)
%   creates a Transfer object, which can hold a transfer function. The idea 
%   is that during a measurement (see Transfer/measure) a stimulus is 
%   played over the DA, and the response is recorded over the AD. Inputs 
%   arguments are
%
%        Q_stim:  quantity of DA signal, e.g., 'Voltage'
%        Q_resp:  quantity of AD signal, e.g., 'Pressure'
%      Ref_stim:  reference for stimulus, e.g. '1 volt RMS'
%      Ref_resp:  reference for response, e.g. '20 uP'
%   dBref_stim:  decibel "unit" for stimulus, e.g. 'dBV'
%   dBref_resp:  decibel "unit" for response, e.g. 'dB SPL'
%
%   The last two inputs, dBref_stim and dBref_resp, are optional. If omitted,
%   they are assigned values 'dB re ...' where ... stands for Ref_stim and
%   Ref_resp, respectively.
%
%   Transfer(), with no input args, returns a void Transfer object. 
%
%   Type "methods transfer" to know what can be done with and to Transfer
%   objects.

if nargin<1, 
    [Q_stim, Q_resp, Ref_stim, Ref_resp, dBref_stim, dBref_resp] = deal([]);
elseif nargin==1 && isstruct(Q_stim), % conversion strct->transfer
    Ttemplate = struct(transfer()); % use template to ensure right field names
    T = structJoin(Ttemplate, Q_stim);
    T = class(T, mfilename);
    return;
end
if nargin<5,
    dBref_stim = ['dB re ' Ref_stim]; 
end
if nargin<6,
    dBref_resp = ['dB re ' Ref_resp]; 
end

% the contents to be filled by methods
Description = '';
CalibParam = [];
Freq = [];
Ztrf = [];
WB_delay_ms = [];
Userdata = struct;

T = CollectInStruct(Q_stim, Q_resp, Ref_stim, Ref_resp, dBref_stim, dBref_resp, ...
    Description, CalibParam, Freq, Ztrf, WB_delay_ms, Userdata);
T = class(T, mfilename);













