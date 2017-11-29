clear all
close all



%alldatain=importdata('C:\Users\Oliver\Downloads\mobile_taxonomy_170410\mobile_taxonomy\MATLABdata_zsolt_new\data\alldata_newest.xlsx');
alldatainMAT=importdata('C:\Users\Oliver\Downloads\mobile_taxonomy_170410\mobile_taxonomy\MATLABdata_Zsolt_Newest\data\datasum.xlsx');
alldatainEFEL = importdata('C:\Users\Oliver\Downloads\Processed\eFEL sum selected.xlsx');
savefileScore = 'C:\Users\Oliver\Desktop\Egyetem\7\Szakdoga\Matlab\Feature_correl.txt';
savefileLabel = 'C:\Users\Oliver\Desktop\Egyetem\7\Szakdoga\Matlab\Feature_corrlabel.txt';
savefileInfo = 'C:\Users\Oliver\Desktop\Egyetem\7\Szakdoga\Matlab\Feature_corrinfo_above_90.txt';
fIDcor = fopen(savefileScore, 'w');
fIDlabel= fopen(savefileLabel, 'w');
fIDinfo= fopen(savefileInfo, 'w');
% celltype=alldatain.textdata(2:end,2);
%experiment=alldatain.textdata(2:end,2);
%numbers = alldatain.data(:,1)';
featureMAT=alldatainMAT.textdata(1,3:end);
alldataMAT=alldatainMAT.data(:,3:end);
numoffeaturesMAT = size(alldataMAT,2);
featureEFEL=alldatainEFEL.textdata(1,3:end);
alldataEFEL=alldatainEFEL.data(:,3:end);
numoffeaturesEFEL = size(alldataEFEL,2);

corrmatrix = cell(numoffeaturesEFEL,numoffeaturesMAT);
corrlabel = cell(numoffeaturesEFEL,numoffeaturesMAT);

for featnumEFEL = 1:numoffeaturesEFEL
    for featnum2MAT = 1:numoffeaturesMAT
        C = corrcoef(alldataEFEL(:,featnumEFEL),alldataMAT(:,featnum2MAT));
        %corrmatrix(featnumEFEL,featnumEFEL) = {C(1,1)};
        corrmatrix(featnumEFEL,featnum2MAT) = {C(1,2)}; %symmetry of the matrix
        %corrmatrix(featnum2MAT,featnum2MAT) = {C(2,2)};
        scoreVal = C(1,2);
        %corrlabel(featnumEFEL,featnumEFEL) = {correl_Label(C(1,1))};
        corrlabel(featnumEFEL,featnum2MAT) = {correl_Label(C(1,2))};
        %corrlabel(featnum2MAT,featnum2MAT) = {correl_Label(C(2,2))};
        
        % Optional: creating files that lists features with a certain correl. label
        scoreLabel = correl_Label(C(1,2));
         if strcmp(scoreLabel,'very strong-') || strcmp(scoreLabel,'very strong+')
             fprintf(fIDinfo, '%s - %s\t%s %.4f \r\n', featureEFEL{featnumEFEL}, featureMAT{featnum2MAT}, scoreLabel, scoreVal);
         end
    end
end
fprintf(fIDcor,'%s\t',' ');
fprintf(fIDlabel,'%s\t',' ');
for ind = 1:size(featureMAT, 2)
    fprintf(fIDcor,'%s\t',featureMAT{ind});
    fprintf(fIDlabel,'%s\t',featureMAT{ind});
end
fprintf(fIDcor,'\r\n');
fprintf(fIDlabel,'\r\n');

for ind1 = 1:size(corrmatrix,1)
    fprintf(fIDcor,'%s\t',featureEFEL{ind1});
    fprintf(fIDlabel,'%s\t',featureEFEL{ind1});
    for ind2 = 1:size(corrmatrix,2)
        fprintf(fIDcor,'%f\t',corrmatrix{ind1,ind2});
        fprintf(fIDlabel,'%s\t', corrlabel{ind1,ind2});
    end
    fprintf(fIDcor,'\r\n');
    fprintf(fIDlabel,'\r\n');
end

fclose(fIDcor);
fclose(fIDlabel);
fclose(fIDinfo);
close all


