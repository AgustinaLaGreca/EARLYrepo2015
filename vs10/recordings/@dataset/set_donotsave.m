function DS = set_donotsave( DS ,donotsave)
% set_donotsave sets the donotsave field is the Datasets ID field

DS = download(DS);
DS.ID.donotsave = donotsave;
DS = upload(DS);
return;

end

