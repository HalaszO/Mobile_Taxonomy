%% MG cluster analysis
[valuenames,fnames]=getvaluenames(datasum,1);
mgvalues={'SAG','Rin','vresting','reobase','aphw', 'apampl','threshold','apwidth','ahpampl','ahpwidth','accomodation','ramp'};
valuenames=mgvalues;
[out,datasum]=getdatafromdatasumr(fnames, datasum, data, iv, valuenames,0);
out.ngfpyrfs=[out.NGFsum;out.PyrL2(2:end,:);out.FSsum(2:end,:)];
clusterstuff=clustermydata(out.ngfpyrfs);

%% distances
tavolsag='seuclidean';
valuenames=AACVSFSVALUES;
[valuenames,fnames]=getvaluenames(datasum,1);
[out,datasum]=getdatafromdatasumr(fnames, datasum, data, iv, valuenames,0);
out.aacmeanvsfs=[out.AACsum(1,:);{'themeanaac'},mat2cell(mean(cell2mat(out.AACsum(2:end,2:end))),1,ones(1,size(out.AACsum,2)-1));out.AACsum(2:end,:);out.FSsum(2:end,:);out.BCsum(2:end,:)];
distances.all=squareform(pdist(cell2mat(out.aacmeanvsfs(2:end,2:end)),tavolsag));%,abs(nanmean(cell2mat(out.aacmeanvsfs(2:end,2:end))))));
valuenames=Sanyivalues;
[out,datasum]=getdatafromdatasumr(fnames, datasum, data, iv, valuenames,0);
out.aacmeanvsfs=[out.AACsum(1,:);{'themeanaac'},mat2cell(mean(cell2mat(out.AACsum(2:end,2:end))),1,ones(1,size(out.AACsum,2)-1));out.AACsum(2:end,:);out.FSsum(2:end,:);out.BCsum(2:end,:)];
distances.sanyi=squareform(pdist(cell2mat(out.aacmeanvsfs(2:end,2:end)),tavolsag));%,abs(nanmean(cell2mat(out.aacmeanvsfs(2:end,2:end))))));
valuenames=Marcivalues;
[out,datasum]=getdatafromdatasumr(fnames, datasum, data, iv, valuenames,0);
out.aacmeanvsfs=[out.AACsum(1,:);{'themeanaac'},mat2cell(mean(cell2mat(out.AACsum(2:end,2:end))),1,ones(1,size(out.AACsum,2)-1));out.AACsum(2:end,:);out.FSsum(2:end,:);out.BCsum(2:end,:)];
distances.marci=squareform(pdist(cell2mat(out.aacmeanvsfs(2:end,2:end)),tavolsag));%,abs(nanmean(cell2mat(out.aacmeanvsfs(2:end,2:end))))));

%% plot distances
toplot={'all','sanyi','marci'};
figure;
for i=1:length(toplot)
    subplot(length(toplot),1,i);
    plot(distances.(char(toplot(i)))(:,1),'o');
    legend(toplot(i));
end
%% get significant values for cluster analysis
ezekbiztosannemkellenek={'reobasesweep','steadysweep','steadyAPnum','burstinterval','bursting'};
%temp=histmysummeddata(out,{'Allsum','AACsum','FSsum','BCsum'});
[valuenames,fnames]=getvaluenames(datasum,1);
[out,datasum]=getdatafromdatasumr(fnames, datasum, data, iv, valuenames,0);
temp=histmysummeddata(out,{'Allsum','AACsum','FSsum','BCsum'});
significantvalues=temp(:,1)';
for i=1:length(ezekbiztosannemkellenek)
    significantvalues(strcmp(significantvalues,char(ezekbiztosannemkellenek(i))))=[];
end
[out,datasum]=getdatafromdatasumr(fnames, datasum, data, iv, significantvalues,0);
out.aacvsfs=[out.AACsum;out.FSsum(2:end,:);out.BCsum(2:end,:)];
clusterstuff=clustermydata(out.aacvsfs);
%% cluster
valuenames=AACVSFSVALUES;
[out,datasum]=getdatafromdatasumr(fnames, datasum, data, iv, valuenames,0);
out.aacvsfs=[out.AACsum;out.FSsum(2:end,:);out.BCsum(2:end,:)];
clusterstuff=clustermydata(out.aacvsfs);

%% principális komponenses fos
aacnum=size(out.AACsum,1)-1;
fsnum=size(out.FSsum,1)-1;
bcnum=size(out.BCsum,1)-1;
[out,datasum]=getdatafromdatasumr(fnames, datasum, data, iv, Marcivalues,0);
out.aacvsfs=[out.AACsum;out.FSsum(2:end,:);out.BCsum(2:end,:)];
[pc,score,latent,tsquare] = princomp(zscore(cell2mat(out.aacvsfs(2:end,2:end))));
figure;
% biplot(pc(:,1:3),'Scores',score(:,1:3));%,'VarLabels',{'X1' 'X2' 'X3' 'X4'});
hold on;
plot3(score(1:aacnum,1),score(1:aacnum,2),score(1:aacnum,3),'ro')
plot3(score(aacnum+1:aacnum+fsnum,1),score(aacnum+1:aacnum+fsnum,2),score(aacnum+1:aacnum+fsnum,3),'go')
plot3(score(aacnum+fsnum+1:end,1),score(aacnum+fsnum+1:end,2),score(aacnum+fsnum+1:end,3),'bo')

%% test values
testvalues=mgvalues;
[valuenames,fnames]=getvaluenames(datasum,1);
goodvalues=[];
for i=1:length(testvalues)
    goodvalues(i)=any(strcmp(valuenames,char(testvalues(i))));
end
if min(goodvalues)==0
    ['SZOPI']
    goodvalues
end