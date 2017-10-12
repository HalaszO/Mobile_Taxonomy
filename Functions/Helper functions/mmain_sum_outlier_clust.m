%% data
valuenamesidontcare={'noiselevel','filterednoiselevel','v0','Rinstd','reobasesweep','steadysweep','steadyAPnum','goodsteadyAPnum','taumin','taunewbest','burstinterval','ramp','rheobaseramp','bursting'};
Marcivalues={'Rin','reobase','vresting','maxapnum','aphw','apampl','threshold','apstartenddiff','dvmax','dvmaxv','dvmin','dvminv','ahpampl','ahpwidth','firstapmaxtime','rheobaserampnew'};
Sanyivalues={'Rin','reobase','aphw','apampl','dvmax','dvminv','ahpamplnew','ahpwidthnew','maxapnum'};
mgvalues={'SAG','Rin','vresting','reobase','aphw', 'apampl','threshold','apwidth','ahpampl','ahpwidth','accomodation','ramp'};
Sanyivalues={'steadyapampl','steadyahpampl','dvmin','dvmax','reobase','Rin','apampl'};
AACVSFSVALUES={'Rin','apampl','dvmax','dvmin','reobase','ahpampl','ahpwidth' };
%preprocessed_folder='/data/Userek/marci/analizis/_neuron taxonomy/FS_vs_AAC/processed_DATA_IV/';
preprocessed_folder='/data/Userek/marci/analizis/_neuron taxonomy/RAT harvested cells/processed_DATA_IV/';
makeimagesofclusters=1;
getallvalues=1;
plotandsaveivs=0;
keepdataandivfiles=1;

%% beolvasás
cd(preprocessed_folder);
fnames=dir;
clear temp;
for i=1:length(fnames)
    temp(i)=fnames(i).isdir;
end
fnames(find(temp))=[];
clear temp;

a=fnames;
clear fnames;
progressbar('reading and analysing files');
for i=1:length(a)
    if any(strfind(a(i).name,'.mat'))
        loadtemp=load([preprocessed_folder,a(i).name]);
        experiment=(char(fieldnames(loadtemp.iv)));
        data.(experiment)=loadtemp.data.(experiment);
        iv.(experiment)=loadtemp.iv.(experiment);
        datasum.(experiment)= calculateelfiz(iv.(experiment),data.(experiment));
        if plotandsaveivs==1
            plotiv(experiment, iv, datasum, data, 1, 0);
            saveas(gcf, [experiment,'.jpg']);
            %close all;
            close(gcf);
        end
        if ~(keepdataandivfiles==1)
            data=rmfield(data,experiment);
            iv=rmfield(iv,experiment);
        end
    end
    progressbar(i/length(a));
end
clear a temp loadtemp;

%% rest
% progressbar('analysing data');
% cellnames=fieldnames(iv);
% for i=1:length(cellnames)
% %     fname=cell2mat(fnames(i));
% %     cellname=fname(1:findstr(fname,'.')-1)
% %     cellname=['elfiz',cellname];
%     cellname=char(cellnames(i))
%
%     progressbar(i/length(cellnames));
% end
%experiments=fieldnames(data);
% for i=1:length(experiments)
%     experiment=char(experiments(i))
%     dlmwrite([experiment,'.txt'],datasum.(experiment).apnum, '\t');
% end
% datasum=getagefromxls(datasum);


[valuenames,fnames]=getvaluenames(datasum,1);
if ~(getallvalues==1)
    %     for i=1:length(valuenamesidontcare)  %Ha nem kell minden valuename %
    %         valuenames(ismember(valuenames,valuenamesidontcare(i)))=[];
    %     end
    valuenames=valuenamesineed;
end
[out,datasum]=getdatafromdatasumr(fnames, datasum, data, iv, valuenames,0);
% [valuenames,fnames]=getvaluenames(datasum,1);

% % outliers=getfakenoutliers(out.Allsum,iv, datasum, data);
% pyroutliers=getfakenoutliers(out.Pyrsum,iv, datasum, data);
% pyroutliers=listarendezo(pyroutliers);
% ngfoutliers=getfakenoutliers(out.NGFsum,iv, datasum, data);
% ngfoutliers=listarendezo(ngfoutliers);
% fsoutliers=getfakenoutliers(out.FSsum,iv, datasum, data);
% fsoutliers=listarendezo(fsoutliers);
%  L2pyroutliers=getfakenoutliers(out.PyrL2,iv, datasum, data);
%  L2pyroutliers=listarendezo(L2pyroutliers);
%  L5pyroutliers=getfakenoutliers(out.PyrL5,iv, datasum, data);
%  L5pyroutliers=listarendezo(L5pyroutliers);
%
% mcorrelatemydata(out);
% for j=1:size(L2pyroutliers,1)
%     fnames(ismember(fnames,L2pyroutliers(j,1)))=[];
% end
% % out=getdatafromdatasumr(fnames, datasum, data, iv, valuenames);
clusterstuff=clustermydata(out.Allsum);
%% sanyi féle plottolás
clusters=fieldnames(clusterstuff);
if makeimagesofclusters==1
    for i=1:length(clusters)
        cluster=char(clusters(i));
        cells=clusterstuff.(cluster);
        for j=1:length(cells)
            cell=char(cells(j));
            plotiv(['elfiz',cell], iv, datasum, data, 0, 2500);
            saveas(gcf, [cluster,'_',cell,'.tiff'],'tiffn');
            close(gcf);cclose all
            clear all
            
%             figure(111);
%             hold all;;
%             plot(datasum.(['elfiz',cell]).reobase:20:20*(length(datasum.(['elfiz',cell]).apnumfromreobase)-1)+datasum.(['elfiz',cell]).reobase,datasum.(['elfiz',cell]).apnumfromreobase,'-o')
        end
%         ylim([0 100])
%         xlim([0 1500])
%         legend(gca,cells);
%         saveas(gcf, [cluster,'.jpg']);
%         close all;
    end
end


