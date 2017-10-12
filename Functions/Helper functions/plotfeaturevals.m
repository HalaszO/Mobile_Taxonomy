clear all
close all

alldatain=importdata('/home/heivi/mobile_taxonomy/MATLABdata/data/allhippo/allmatdatanew.csv');

celltype=alldatain.textdata(2:end,1);
experiment=alldatain.textdata(2:end,2);
feature=alldatain.textdata(1,3:end);

celltypes=unique(celltype);
include=celltypes(4:8);

for fnum=1:length(feature),
    for ct=1:length(include),
        cellids=strmatch(include(ct),celltype,'exact');
        data=alldatain.data(cellids,fnum);
        plot(ct,data,'kx')
        hold on
        ylabel(feature(fnum))
        ax=axis;
        axis([0 length(include)+1 ax(3) ax(4)])
    end
    hold off
    pause
end
