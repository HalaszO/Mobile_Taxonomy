
%alldatain=importdata('C:\Users\Oliver\Desktop\Egyetem\7\Szakdoga\Matlab\eFEL-MATLAB correlation\eFEL MATLAB correlation.xlsx');
alldatain = importdata('C:\Users\Oliver\Desktop\Egyetem\7\Szakdoga\Matlab\eFEL correlation\efel cor.xlsx');
%alldatain = importdata('C:\Users\Oliver\Desktop\Egyetem\7\Szakdoga\Matlab\MATLAB correlation\MATLAB correl.xlsx');
%experiment=alldatain.textdata(2:end,4);
MATLABfeatures = alldatain.textdata.Sheet1(1,2:end);
eFELfeatures = alldatain.textdata.Sheet1(2:end,1);
data = alldatain.data.Sheet1(2:end,2:end);
picsavefolder = 'C:\Users\Oliver\Desktop\Egyetem\7\Szakdoga\Correlation\';

% weakfind = (data >= -0.35 & data <= 0.35);
% weak = sum(sum(weakfind));
% all = 74*61;
% percent = weak/all*100

figure(1)
pcolor(data)
colormap jet
colorbar
title('Correlation')
xlabel('MATLAB features')
ylabel('eFEL features')
%saveas(gcf,[picsavefolder,'MATLAB self correlation.png']);
%saveas(gcf,[picsavefolder,'MATLAB eFEL correlation.fig']);