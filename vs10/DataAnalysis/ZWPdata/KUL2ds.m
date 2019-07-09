function d = KUL2ds(FN)
% KUL2ds - converter KUL => Rdam EARLY compatible dataset
%    KUL2ds('file') reads file and writes file_rdam

if ischar(FN),
    FN = strrep(FN, '.mat', '');
    load(FN); % loads var named d
else,
    d = FN;
end
d = struct(d);
for ii=1:numel(d.Rec.RecordInstr),
    g = d.Rec.RecordInstr(ii).grabber;
    g = fhandle(func2str(g));
    d.Rec.RecordInstr(ii).grabber = g;
end

if isstruct(d.Stim.Presentation),
    d.Stim.Presentation = rmfield(d.Stim.Presentation, 'ReducedStorage');
end

d.Data.RX6_digital_1 = struct(d.Data.RX6_digital_1);
d.Data.RX6_digital_1.grabbeddata = struct(d.Data.RX6_digital_1.grabbeddata);
d.Data.RX6_digital_1.grabbeddata.Supplier = struct(d.Data.RX6_digital_1.grabbeddata.Supplier);
d.Data.RX6_digital_1.grabbeddata.Supplier.Hardware = struct(d.Data.RX6_digital_1.grabbeddata.Supplier.Hardware);
d.Data.RX6_digital_1.grabbeddata.Supplier.Hardware.ResetDAC_trigger = [];
d.Data.RX6_digital_1.grabbeddata.Supplier.Hardware.CalibStamping_trigger = [];
d.Data.RX6_digital_1.grabbeddata.Supplier.Hardware.ResetStamping_trigger = [];
d.Data.RX6_digital_1.grabbeddata.Supplier.Hardware.CalibADC_trigger = [];
d.Data.RX6_digital_1.grabbeddata.Supplier.Hardware.ResetADC_trigger = [];
d.Data.RX6_digital_1.grabbeddata.Supplier.datagrabber = struct(d.Data.RX6_digital_1.grabbeddata.Supplier.datagrabber);
d.Data.RX6_digital_1.grabbeddata.Supplier.datagrabber.RecInstr.grabber = fhandle(func2str(d.Data.RX6_digital_1.grabbeddata.Supplier.datagrabber.RecInstr.grabber));
d.Data.RX6_digital_1.grabbeddata.Supplier.datagrabber.DS = [];
d.Data.RX6_digital_1.grabbeddata.Supplier.datagrabber.SourceDevice.ID.RecInstr.grabber= [];
%d.Data.RX6_digital_1.grabbeddata.Supplier.datagrabber = [];

% d.Data.RX6_analog_1= struct(d.Data.RX6_analog_1);
d.Data = rmfield(d.Data, 'RX6_analog_1');

% undo the vorseid obj=>struct conversions
d.Data.RX6_digital_1.grabbeddata.Supplier.datagrabber = struct2obj(d.Data.RX6_digital_1.grabbeddata.Supplier.datagrabber, 'datagrabber');
d.Data.RX6_digital_1.grabbeddata.Supplier.Hardware = struct2obj(d.Data.RX6_digital_1.grabbeddata.Supplier.Hardware, 'rechardware');
d.Data.RX6_digital_1.grabbeddata.Supplier = struct2obj(d.Data.RX6_digital_1.grabbeddata.Supplier, 'grabevents');
d.Data.RX6_digital_1.grabbeddata = struct2obj(d.Data.RX6_digital_1.grabbeddata, 'grabbeddata');
d.Data.RX6_digital_1 = struct2obj(d.Data.RX6_digital_1, 'eventdata');
d.Stim.Presentation = struct2obj(d.Stim.Presentation, 'stimpresentx');
d = struct2obj(d, 'dataset');

if ischar(FN),
    save([FN '_rdam'], 'd', '-v6');
end















