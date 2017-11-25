clear all
close all

%%%%alldatain = importdata('C:\Users\Oliver\Desktop\Egyetem\7\Szakdoga\HBP_processed_all_grouped_new\data\eFEL processed HBP grouped.xlsx');
alldatain = importdata('C:\Users\Oliver\Desktop\Egyetem\7\Szakdoga\HBP_processed_all_grouped_new\data\MATLAB HBP all processed.xlsx');
% celltype=alldatain.textdata(2:end,2);
experiment=alldatain.textdata(2:end,5);
%expcount = alldatain.data(1:end,1)';
expcount = alldatain.data(1:77,1)';
features=alldatain.textdata(1,5:end)';
% alldata=alldatain.data(1:77,5:end);
 %alldata=alldatain.data(1:133,5:end);
alldata=alldatain.data(:,5:end);


picsavefolder = ('C:\Users\Oliver\Desktop\Egyetem\7\Szakdoga\HBP_pics\HBP_data_all_grouped\MATLAB results\ALL');
%picsavefolder = ('C:\Users\Oliver\Desktop\Egyetem\7\Szakdoga\HBP_pics\HBP_data_all_grouped\MATLAB results\PC_CCK_PRINC');
if ~exist(picsavefolder,'dir') 
    mkdir(picsavefolder); 
end

%%%%%%%%%%%%%%%
% PCA
%%%%%%%%%%%%%%%

stdr = std(alldata);
alldata_norm = alldata./repmat(stdr,size(alldata,1),1);
[coefs,scores,variances,t2] = princomp(alldata_norm);

%%%%%%%%%%%
% Saving coefs
%%%%%%%%%%

file_coef = fopen([picsavefolder,'\MATLAB_coefs_1.txt'],'w');
for ind = 1:size(coefs, 1)
    fprintf(file_coef,'%s\t%f\r\n',features{ind},coefs(ind,1));
end
fclose(file_coef);

file_coef = fopen([picsavefolder,'\MATLAB_coefs_2.txt'],'w');
for ind = 1:size(coefs, 1)
    fprintf(file_coef,'%s\t%f\r\n',features{ind},coefs(ind,2));
end
fclose(file_coef);

%%%%%%%%%%%%%%%
% Plot PC-2 scores
%%%%%%%%%%%%%

figure
hold on
plot(scores(1:57,1),scores(1:57,2),'x', 'MarkerSize',10, 'MarkerEdgeColor','blue')
plot(scores(58:77,1),scores(58:77,2),'x','MarkerSize',10, 'MarkerEdgeColor',[1, 0.6, 0])
plot(scores(78:133,1),scores(78:133,2),'x','MarkerSize',10, 'MarkerEdgeColor',[0.15, 1, 0])
plot(scores(134:346,1),scores(134:346,2),'x','MarkerSize',6, 'MarkerEdgeColor',[0.73 0.73 0.73])
xlabel('First PC')
ylabel('Second PC')
title('Data in first two PC space')
legend('CCK+', 'PV+', 'Principal cell', 'Other and unknown')
saveas(gcf,[picsavefolder,'\second dotplot.png']);
hold off

%%%%%%%%%%%%%%%%%
% PC variance
%%%%%%%%%%
percent_explained = 100*variances/sum(variances);
figure
pareto(percent_explained)
title('Variance of the PCs')
xlabel('Principal component')
ylabel('Variance')
saveas(gcf,[picsavefolder,'\HBP principal components var.png']);

%%%%%%%%%%%%%%%%%
% Tree construction
%%%%%%%%%%%%%%%
%2
scores2=scores(:,1:2);
dist2=pdist(scores2,'euclidean');   % euclidean dist
link2=linkage(dist2,'ward'); % hierarchical cluster tree, using Inner squared distance
figure
dendrogram(link2, 0, 'colorThreshold','default')
title('Dendogram using 2 PCs')
xlabel('Data points')
ylabel('Distance')
set(gca,'XTickMode','auto','XTickLabelMode','auto');
saveas(gcf,[picsavefolder,'\HBP dendo coefs 2.png']);

%4
scores4=scores(:,1:4);
dist4=pdist(scores4,'euclidean');
link4=linkage(dist4,'ward');
figure
dendrogram(link4,0, 'colorThreshold','default')
title('Dendogram using 4 PCs')
xlabel('Data points')
ylabel('Distance')
set(gca,'XTickMode','auto','XTickLabelMode','auto');
saveas(gcf,[picsavefolder,'\HBP dendo coefs 4.png']);

%6
scores6=scores(:,1:6);
dist6=pdist(scores6,'euclidean');
link6=linkage(dist6,'ward');
figure
dendrogram(link6,0, 'colorThreshold','default')
title('Dendogram using 6 PCs')
xlabel('Data points')
ylabel('Distance')
set(gca,'XTickMode','auto','XTickLabelMode','auto');
saveas(gcf,[picsavefolder,'\HBP dendo coefs 6.png']);

%8
scores8=scores(:,1:8);
dist8=pdist(scores8,'euclidean');
link8=linkage(dist8,'ward');
figure
dendrogram(link8,0, 'colorThreshold','default')
title('Dendogram using 8 PCs')
xlabel('Data points')
ylabel('Distance')
set(gca,'XTickMode','auto','XTickLabelMode','auto');
saveas(gcf,[picsavefolder,'\HBP dendo coefs 8.png']);

ccc = zeros(1,4);
ccc(1) = cophenet(link2, dist2);
ccc(2) = cophenet(link4, dist4);
ccc(3) = cophenet(link6, dist6);
ccc(4) = cophenet(link8, dist8);

file_coef = fopen([picsavefolder,'\coph_coefs.txt'],'w');
fprintf(file_coef,'%f\t',ccc);
fclose(file_coef);

%%%%%%%%%%%%%%%%%%%%%
% Inconsistency measurement
%%%%%%%%%%

in2 = inconsistent(link2);
in4 = inconsistent(link4);
in6 = inconsistent(link6);
in8 = inconsistent(link8);
% maxin8 = max(in4(:,4));
maxin6 = max(in6(:,4));
maxin4 = max(in8(:,4));
% maxin2 = max(in2(:,4));

%%%%%%%%%%%%%%%%%%%%%%%%
% Agglomerative clusters from inconsistency
%%%%%%%%%%%%
figure
c = cluster(link6,'cutoff',maxin6-eps);
scatter(scores2(:,1),scores2(:,2),30,c,'filled')
title('Clustering based on link inconsistency')
xlabel('First PC')
ylabel('Second PC')
saveas(gcf,[picsavefolder,'\HBP cluster inconsist coefs 8.png']);

figure
silhouette(scores6,c,'Euclidean')
title('Silhouette values of inconsistency clustering')
xlabel('Silhouette values')
ylabel('Clusters')
saveas(gcf,[picsavefolder,'\HBP silhouette inconsist 8.png']);


%%%%%%%%%%%%%%%%%%%%
% Evaluation of clustering
%%%%%%%%%

% eva2 = evalclusters(link2,'linkage','CalinskiHarabasz')
% eva4 = evalclusters(link4,'linkage','CalinskiHarabasz')
% eva6 = evalclusters(link6,'linkage','')
eva4 = evalclusters(scores6,'linkage','Silhouette', 'KList',[1:10])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Agglomerative clusters with optimal cluster num
%%%%%%%%%%%%
figure
c = cluster(link6,'maxclust',eva4.OptimalK);
scatter(scores2(:,1),scores2(:,2),30,c,'filled')
title('Clustering based on optimal cluster number')
xlabel('First PC')
ylabel('Second PC')
saveas(gcf,[picsavefolder,'\HBP cluster optim num coefs 8.png']);

figure
silhouette(scores6,c,'Euclidean')
title('Silhouette values of optimal clustering')
xlabel('Silhouette values')
ylabel('Clusters')
saveas(gcf,[picsavefolder,'\HBP silhouette optim 8.png']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fuzzy c-means clustering
%%%%%%%%%%%%%

[~,U] = fcm(scores6,2);
%Classify each data point into the cluster with the largest membership value.
maxU = max(U);
index1 = find(U(1,:) == maxU);
index2 = find(U(2,:) == maxU);

fuz_sil = zeros(size(alldata, 1),1);
fuz_sil(index1) = 1;
fuz_sil(index2) = 2;

%Plotting
figure
hold on
plot(scores2(index1,1),scores2(index1,2),'+b')
plot(scores2(index2,1),scores2(index2,2),'+r')
xlabel('First PC')
ylabel('Second PC')
title('Fuzzy-clustering in the first two PC space')
%legend('Cluster 1', 'Cluster 2')
saveas(gcf,[picsavefolder, '\HBP fuzzy cluster 2.png']);
hold off

figure

ss2 = silhouette(scores6,fuz_sil,'Euclidean');
silhouette(scores6,fuz_sil,'Euclidean')
title('Silhouette values of fuzzy clustering')
xlabel('Silhouette values')
ylabel('Clusters')
saveas(gcf,[picsavefolder,'\HBP silhouette fuzzy 2.png']);

ss2mean = mean(ss2);
ss2std = std(ss2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fuzzy c-means clustering
%%%%%%%%%%%%%

[~,U] = fcm(scores6,3);
%Classify each data point into the cluster with the largest membership value.
maxU = max(U);
index1 = find(U(1,:) == maxU);
index2 = find(U(2,:) == maxU);
index3 = find(U(3,:) == maxU);


fuz_sil = zeros(size(alldata, 1),1);
fuz_sil(index1) = 1;
fuz_sil(index2) = 2;
fuz_sil(index3) = 3;

%Plotting
figure
hold on
plot(scores2(index1,1),scores2(index1,2),'+b')
plot(scores2(index2,1),scores2(index2,2),'+r')
plot(scores2(index3,1),scores2(index3,2),'+g')
xlabel('First PC')
ylabel('Second PC')
title('Fuzzy-clustering in the first two PC space')
%legend('Cluster 1', 'Cluster 2')
saveas(gcf,[picsavefolder, '\HBP fuzzy cluster 3.png']);
hold off

figure
ss2 = silhouette(scores6,fuz_sil,'Euclidean');
silhouette(scores6,fuz_sil,'Euclidean')
title('Silhouette values of fuzzy clustering')
xlabel('Silhouette values')
ylabel('Clusters')
saveas(gcf,[picsavefolder,'\HBP silhouette fuzzy 3.png']);

ss3mean = mean(ss2);
ss3std = std(ss2);

% incor_pv = sum(index2 < 58);
% incor_cck = sum(index1 > 57);
% correct_clust = 1 - ((incor_pv+incor_cck)/size(expcount,2));
 
% incor_princ = sum(index2 < 78);
% incor_int = sum(index1 > 77);
% correct_clust = 1 - ((incor_princ+incor_int)/size(expcount,2));
