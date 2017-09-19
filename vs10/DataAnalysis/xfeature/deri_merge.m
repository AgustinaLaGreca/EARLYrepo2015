function EV = deri_merge(AP,EPSP,D,DS)
% deri_merge - Unify and merge AP and EPSP fields of the deristat results
%   This function takes the deristat results of one condition of one recording
%   of one experiment and merges them into a unified format. The advantage
%   is, that now all events are contained, regardless of their type
%   (i.e. successful or failed event).
%
%   syntax: EV = deri_merge(AP,EPSP,D,DS)
%
%   A number of new fields per event are introduced. They will be listed
%   here:
%   EV.isAP... logical that indicates, whether event is AP or not
%   EV.isinStim... logical that indicates, whether event is during stim
%   EV.StimNum...which stimulus was it (0 == none)
%   EV.priorHistory...how many events of the same type happened before in
%                       that stimulus
%   EV.ExpID... SGSR ExpID
%   EV.RecID... SGSR RecID
%   EV.Condition... Condition in that recording
%   EV.xval... xval of that condition
%   EV.xunit... xval of that condition
%   EV.yval... yval of that recording
%   EV.yunit... yunit of that recording
%   EV.carfreq... Carrier Frequency
%   EV.repdur... duration of one repetition
%   EV.burstdur... duration of one stimulus presentation
%   EV.stimtype... type of stimulus (FS or SPL, etc...)
%
%   See also readTKabf, deristats

%%%%%%%%%%%%some preparation
EV = struct;
NAP = numel(AP);
NEPSP = numel(EPSP);
Nall = NAP + NEPSP;

%%%%%%%%%%%%Collect stimulus information
if strncmpi(DS.Stimtype,'FS',2),
    yunit = 'db SPL';
    yval = DS.SPL(1);
    carfreq = DS.carfreq(D.icond);
    Nrep = DS.param.stimcntrl.repcount;
elseif strcmp(DS.Stimtype,'SPL'),
    yunit = 'Hz';
    yval = DS.fcar;
    carfreq = DS.carfreq;
    Nrep = DS.param.stimcntrl.repcount;
elseif strcmp(DS.Stimtype,'TTS'),
    yunit = 'Hz(suppr)';
    yval = DS.param.suppfreq;
    carfreq = DS.carfreq;
    Nrep = DS.reps;
end;
repdur = D.sweepDur_ms;
stimdur = DS.burstdur; 
stimtype = DS.StimType;


%%%%%%%%%%%%%%%deal result fields
%divider___
%isAP?
[EV(1:NAP).isAP] = deal(true);%APs are APs
[EV(NAP+1:Nall).isAP] = deal(false);%EPSPs are no APs
[EV.a____c_o_m_m_o_n____m_e_t_r_i_c_s____] = deal('c_o_m_m_o_n____m_e_t_r_i_c_s__');
%istart & iend (Index of first and last sample that belongs to the event)
[EV(1:NAP).istart] = deal(AP.istart);
[EV(NAP+1:Nall).istart] = deal(EPSP.istart);
[EV(1:NAP).iend] = deal(AP.iend);
[EV(NAP+1:Nall).iend] = deal(EPSP.iend);
%t0 & t1 (Time of event start and end)
[EV(1:NAP).t0] = deal(AP.t0);
[EV(NAP+1:Nall).t0] = deal(EPSP.t0);
[EV(1:NAP).t1] = deal(AP.t1);
[EV(NAP+1:Nall).t1] = deal(EPSP.t1);
%event total duration in ms (EPSPs and APs have that)
[EV(1:NAP).dur] = DealElements([AP.t1] - [AP.t0]);
[EV(NAP+1:Nall).dur] = deal(EPSP.dur);
%Intervals to previous spike and event (EPSPs and APs have that)
[EV(1:NAP).itvAP] = deal(AP.itvAP);
[EV(NAP+1:Nall).itvAP] = deal(EPSP.itvAP);
[EV(1:NAP).itvEvent] = deal(AP.itvEvent);
[EV(NAP+1:Nall).itvEvent] = deal(EPSP.itvEvent); 
%divider___
[EV.b______A_P__m_e_t_r_i_c_s____] = deal('A_P__m_e_t_r_i_c_s____');
%Spike peak timing and voltage (only APs have that)
[EV(1:NAP).Vpeak] = deal(AP.Vpeak);
[EV(NAP+1:Nall).Vpeak] = deal(NaN);
[EV(1:NAP).tpeak] = deal(AP.tpeak);
[EV(NAP+1:Nall).tpeak] = deal(NaN);
%Spike downslope maximum dV and timing (only APs have that)
[EV(1:NAP).maxDownSlope] = deal(AP.maxDownSlope);
[EV(NAP+1:Nall).maxDownSlope] = deal(NaN);
[EV(1:NAP).tmaxDownSlope] = deal(AP.tmaxDownSlope);
[EV(NAP+1:Nall).tmaxDownSlope] = deal(NaN);
%Spike dip maximum V and timing (only APs have that)
[EV(1:NAP).Vdip] = deal(AP.Vdip);
[EV(NAP+1:Nall).Vdip] = deal(NaN);
[EV(1:NAP).tdip] = deal(AP.tdip);
[EV(NAP+1:Nall).tdip] = deal(NaN);
%Spike PostUpslope maximum dV and timing (only APs have that)
[EV(1:NAP).maxPostUpslope] = deal(AP.maxPostUpslope);
[EV(NAP+1:Nall).maxPostUpslope] = deal(NaN);
[EV(1:NAP).tmaxPostUpslope] = deal(AP.tmaxPostUpslope);
[EV(NAP+1:Nall).tmaxPostUpslope] = deal(NaN);
%EPSP-Spike inflexion timing (only APs have that)
[EV(1:NAP).tInflex] = deal(AP.tInflex);
[EV(NAP+1:Nall).tInflex] = deal(NaN);
%Spike monoUpslope logical (only APs have that)
[EV(1:NAP).monoUpslope] = deal(AP.monoUpslope);
[EV(NAP+1:Nall).monoUpslope] = deal(NaN);
%Spike PreUpslope maximum dV and timing (only APs have that)
[EV(1:NAP).maxPreUpslope] = deal(AP.maxPreUpslope);
[EV(NAP+1:Nall).maxPreUpslope] = deal(NaN);
[EV(1:NAP).tmaxPreUpslope] = deal(AP.tmaxPreUpslope);
[EV(NAP+1:Nall).tmaxPreUpslope] = deal(NaN);
%divider___
[EV.c______E_P_S_P__m_e_t_r_i_c_s____] = deal('E_P_S_P__m_e_t_r_i_c_s____');
%event start and end voltage (only EPSPs have that)
[EV(1:NAP).y0] = deal(NaN);
[EV(1:NAP).y1] = deal(NaN);
[EV(NAP+1:Nall).y0] = deal(EPSP.y0);
[EV(NAP+1:Nall).y1] = deal(EPSP.y1);
%EPSP peak timing and voltage (EPSPs and APs have that)
[EV(1:NAP).VpeakEPSP] = deal(AP.VpeakEPSP);
[EV(NAP+1:Nall).VpeakEPSP] = deal(EPSP.Vpeak);
[EV(1:NAP).tpeakEPSP] = deal(AP.tpeakEPSP);
[EV(NAP+1:Nall).tpeakEPSP] = deal(EPSP.tpeak);
%EPSP maxRate timing and voltage (EPSPs and APs have that)
[EV(1:NAP).maxRateEPSP] = deal(AP.EPSPmaxRate);
[EV(NAP+1:Nall).maxRateEPSP] = deal(EPSP.maxUpRate);
[EV(1:NAP).tmaxRateEPSP] = deal(AP.tEPSPmaxRate);
[EV(NAP+1:Nall).tmaxRateEPSP] = deal(EPSP.tmaxUpRate);
%EPSP dip voltage and timing (only EPSPs have that)
[EV(1:NAP).VdipEPSP] = deal(NaN);
[EV(NAP+1:Nall).VdipEPSP] = deal(EPSP.Vdip);
[EV(1:NAP).tdipEPSP] = deal(NaN);
[EV(NAP+1:Nall).tdipEPSP] = deal(EPSP.tdip);
%EPSP downslope max dV and timing (only EPSPs have that)
[EV(1:NAP).maxDownRateEPSP] = deal(NaN);
[EV(NAP+1:Nall).maxDownRateEPSP] = deal(EPSP.maxDownRate);
[EV(1:NAP).tmaxDownRateEPSP] = deal(NaN);
[EV(NAP+1:Nall).tmaxDownRateEPSP] = deal(EPSP.tmaxDownRate);
%event is embedded in a APpostUpflank logical? (only EPSPs have that)
[EV(1:NAP).isAPpostUpFlank] = deal(NaN);
[EV(NAP+1:Nall).isAPpostUpFlank] = deal(EPSP.isAPpostUpFlank);
%%%%%%%%%%%%%%%deal bookkeeping fields
%divider___
[EV.d______e_v_e_n_t__h_i_s_t_o_r_y____] = deal('e_v_e_n_t__h_i_s_t_o_r_y____');
%context and "local history" of the event
[a,b,c] = LOCAL_isinstim(repdur,stimdur,Nrep,[AP.tpeak]);
if NAP > 0,
[EV(1:NAP).isinStim] = DealElements(a);%AP t0 in stim?
[EV(1:NAP).StimNum] = DealElements(b);%which stim?
[EV(1:NAP).priorHistory] = DealElements(c);%how many APs happened before?
end;
[a,b,c] = LOCAL_isinstim(repdur,stimdur,Nrep,[EPSP.tpeak]);
if NEPSP > 0,
[EV(NAP+1:Nall).isinStim] = DealElements(a);%EPSP t0 in stim?
[EV(NAP+1:Nall).StimNum] = DealElements(b);%which stim?
[EV(NAP+1:Nall).priorHistory] = DealElements(c);%how many EPSPs happened before?
end;
%divider___
[EV.e______s_t_i_m_u_l_u_s__i_n_f_o____] = deal('s_t_i_m_u_l_u_s__i_n_f_o____');
%SGSR IDs...
[EV(1:NAP).ExpID] = deal(D.ExpID);%SGSR Experiment ID
[EV(NAP+1:Nall).ExpID] = deal(D.ExpID);%SGSR Experiment ID
[EV(1:NAP).RecID] = deal(D.RecID);%SGSR Recording ID
[EV(NAP+1:Nall).RecID] = deal(D.RecID);%SGSR Recording ID
[EV(1:NAP).Condition] = deal(D.icond);%Number of Condition within the recording
[EV(NAP+1:Nall).Condition] = deal(D.icond);%Number of Condition within the recording
%stimulus and sweep information
[EV(1:NAP).repdur] = deal(repdur);%length of one repetition (corrected) in ms
[EV(NAP+1:Nall).repdur] = deal(repdur);%length of one repetition (corrected) in ms
[EV(1:NAP).stimtype] = deal(stimtype);%SGSR type of stimulus
[EV(NAP+1:Nall).stimtype] = deal(stimtype);%SGSR type of stimulus
[EV(1:NAP).stimdur] = deal(stimdur);%length of one stimulus presentation in ms
[EV(NAP+1:Nall).stimdur] = deal(stimdur);%length of one stimulus presentation in ms
[EV(1:NAP).xval] = deal(DS.xval(D.icond));%value of Stimulus Parameter, that changes with conditions
[EV(NAP+1:Nall).xval] = deal(DS.xval(D.icond));%value of Stimulus Parameter, that changes with conditions
[EV(1:NAP).xunit] = deal(DS.xunit);%Unit of Stimulus Parameter, that changes with conditions
[EV(NAP+1:Nall).xunit] = deal(DS.xunit);%%Unit of Stimulus Parameter, that changes with conditions
[EV(1:NAP).yval] = deal(yval);%value of constant Stimulus Parameter
[EV(NAP+1:Nall).yval] = deal(yval);%value of constant Stimulus Parameter
[EV(1:NAP).yunit] = deal(yunit);%unit of constant Stimulus Parameter
[EV(NAP+1:Nall).yunit] = deal(yunit);%unit of constant Stimulus Parameter
[EV(1:NAP).carfreq] = deal(carfreq);%carrier frequency (relevant for TTS)
[EV(NAP+1:Nall).carfreq] = deal(carfreq);%carrier frequency (relevant for TTS)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%LOCAL functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [during_presentation,which_stim,history] = LOCAL_isinstim(repdur,stimdur,Nrep,tEv)
%LOCAL_isinstim - find stimulus context for events

%calculate times and set up some variables
onset = repdur:repdur:Nrep*repdur;
offset = onset + stimdur;
NtEv = numel(tEv);
during_presentation = false;
which_stim = 0;
history = 0;

%find stuff
if isempty(NtEv) || NtEv == 0,
    return;
elseif NtEv == 1,%we check only one event
    if sum(betwixt(tEv,onset,offset))~=0,%the spiketime is during one stimulus
        during_presentation = true;
        which_stim = find((betwixt(tEv,onset,offset)) ~= 0);
    end;
elseif NtEv > 1,%we check a vector of events
    for nnn = 1:NtEv,
        if sum(betwixt(tEv(nnn),onset,offset))~=0,
            during_presentation(nnn) = true;
            which_stim(nnn) = find((betwixt(tEv(nnn),onset,offset)) ~= 0);
            history(nnn) = numel(find(which_stim == which_stim(nnn)))-1;
        else
            during_presentation(nnn) = false;
            which_stim(nnn) = 0;
            history(nnn) = 0;
        end
    end;
end
