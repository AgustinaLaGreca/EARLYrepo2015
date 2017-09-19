function St=supplierstatus(GD);
% grabbeddata/supplierstatus - status of supplier at time of getdata.
%    supplierstatus(GD) returns the status GD's supplier (the DataGrabber
%    object supplying data for GD) at the time of the Getdata call that
%    created GD.
%
%    See also Grabevents, grabevents/getdata, Datagrabber.
St = GD.SupStatus;






