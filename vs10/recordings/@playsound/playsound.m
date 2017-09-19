function actP=playsound(Wav, Experiment, Rec, Atten, StimPres, varargin);
% playsound - construct playsound object
%    actP=playsound(Wav, EXP, Atten, SP, Ax) creates a playsound object actP 
%    whose task it is to play auditory stimuli. Inputs:
%        W: waveform object or its address, must be uploaded (pointer is used, no copy)
%       EXP: Experiment definition containing generic stimulus settings
%      Rec: instructions for loading the circuit 
%    Atten: instructions for attenuator settings
%       SP: stimulus presentation object
%       Ax: upload destination for ActP (see Dynamic).
% 
%    Obviously, these different pieces of information must be mutually
%    compatible. They are typically created by stimulus makers like
%    tonestim.
%
%    playsound(W, EXP, Atten, SP, figh, 'Foo') uses upload destination
%    Ax==GUIaccess(figh, 'Foo').
%
%    The default timer interval of playsound objects is 44 ms. Note that
%    this does not affect the playing itself, because that is not handled
%    via the timer. It only affects the time it takes for actP to find out
%    that it is ready playing.
%
%    Type 'methods playsound -full' to see what can be done with actP.
%
%  See also tonestim, sortConditions, action/isready, action/start.

isvoid = 0;
status = 'initialized';
if nargin<1,
    [Wav, Experiment, Rec, Atten, StimPres, DT] = deal([]);
    isvoid = 1;
    status = 'void';
end

DT = 44; % ms default timer interval
if isa(Wav, 'Waveform'), Wav=address(Wav); end
actP = CollectInStruct(Wav, Experiment, Rec, Atten, StimPres, DT);

actP = class(actP, mfilename, action(status, DT));
if ~isvoid, actP=upload(actP,varargin{:}); end









