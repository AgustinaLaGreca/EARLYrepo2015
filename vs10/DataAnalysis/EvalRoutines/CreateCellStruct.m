function CellStruct = CreateCellStruct(ds)
CellStruct = [];
CellNr = ds.ID.iCell;
TestNr = ds.ID.iRecOfCell;
dsID = [num2str(CellNr) '-' num2str(TestNr) '-' ds.Stim.StimType];
CellStruct = CollectInStruct(CellNr, TestNr, dsID);
