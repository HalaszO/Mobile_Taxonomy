clear all
close all

%%% Edit paths! %%%

alldatain = importdata('C:\Users\Oliver\Desktop\Egyetem\7\Szakdoga\Elphys_all_fixed_data_again_2\data\Matlabdata_03_28.xlsx');

features=alldatain.textdata(1,3:end)';
alldata=alldatain.data;
featurecount = size(features, 1);
numbering = 1:size(alldata, 1);


picsavefolder = ('C:\Users\Oliver\Desktop\Egyetem\7\Szakdoga\HBP_pics\features_newest_unlabeled');
if ~exist(picsavefolder,'dir') 
    mkdir(picsavefolder); 
end

%%%%%%%%%%%%%%%
% Plotting
%%%%%%%%%%%%%

for idx = 1:featurecount
    figure
    hold on
    plot(1,alldata(1:49,idx),'x', 'MarkerSize',10, 'MarkerEdgeColor','blue')
    %text(ones(1,49),alldata(1:49,idx),num2cell(numbering(1:49)),'FontSize',9)
    plot(2,alldata(50:71,idx),'x','MarkerSize',10, 'MarkerEdgeColor',[1, 0.6, 0])
    %text(2*ones(1,71-50+1),alldata(50:71,idx),num2cell(numbering(50:71)),'FontSize',9)
    plot(3,alldata(72:127,idx),'x','MarkerSize',10, 'MarkerEdgeColor',[0.15, 1, 0])
    %text(3*ones(1,127-72+1),alldata(72:127,idx),num2cell(numbering(72:127)),'FontSize',9)
    plot(4,alldata(128:367,idx),'x','MarkerSize',6, 'MarkerEdgeColor',[0.73 0.73 0.73])
    %text(4*ones(1,367-128+1),alldata(128:367,idx),num2cell(numbering(128:end)),'FontSize',9)
    xlim([0 5])
    ylim('auto')
    xlabel('Cell types: CCK, PV, Principal, Other')
    ylabel('Feature value')
    title(features{idx})
    %legend('CCK+', 'PV+', 'Principal', 'Other and unknown')
    saveas(gcf,[picsavefolder,'\',features{idx},'.png']);
    hold off
end

