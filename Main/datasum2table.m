clear all

%% data
% valuenamesidontcare={'noiselevel','filterednoiselevel','v0','Rinstd','reobasesweep','steadysweep','steadyAPnum','goodsteadyAPnum','taumin','taunewbest','burstinterval','ramp','rheobaseramp','bursting'};
% Marcivalues={'Rin','reobase','vresting','maxapnum','aphw','apampl','threshold','apstartenddiff','dvmax','dvmaxv','dvmin','dvminv','ahpampl','ahpwidth','firstapmaxtime','rheobaserampnew'};
% Sanyivalues={'Rin','reobase','aphw','apampl','dvmax','dvminv','ahpamplnew','ahpwidthnew','maxapnum'};
% mgvalues={'SAG','Rin','vresting','reobase','aphw', 'apampl','threshold','apwidth','ahpampl','ahpwidth','accomodation','ramp'};
% Sanyivalues={'steadyapampl','steadyahpampl','dvmin','dvmax','reobase','Rin','apampl'};
% AACVSFSVALUES={'Rin','apampl','dvmax','dvmin','reobase','ahpampl','ahpwidth' };
%preprocessed_folder='/data/Userek/marci/analizis/_neuron taxonomy/FS_vs_AAC/processed_DATA_IV/';
%preprocessed_folder='C:\Users\Vivien\Downloads\mobile_taxonomy\MATLABdata_HBP/data/HBP_long/';
preprocessed_folder = uigetdir('','Location of processed data files');
makeimagesofclusters=1;
getallvalues=1;
plotandsaveivs=0;
keepdataandivfiles=1;

%% Reading files
cd(preprocessed_folder);
fnames=dir;
clear temp;
for i=1:length(fnames)
    temp(i)=fnames(i).isdir;
end
fnames(temp)=[];
clear temp;

a=fnames;
clear fnames;
progressbar('Reading and analysing files');

for i=1:length(a)
    if any(strfind(a(i).name,'.mat'))
        loadtemp=load([preprocessed_folder,filesep a(i).name]);
        ivname=(char(fieldnames(loadtemp.iv)));
        cellnames=fieldnames(loadtemp.iv.(ivname));
        for ii=1:length(cellnames)
            cellname=char(cellnames(ii));
            experiment=[ivname,'_',cellname];
            if isfield(loadtemp.data.(ivname),cellname)
                data.(experiment)=loadtemp.data.(ivname).(cellname);
                iv.(experiment)=loadtemp.iv.(ivname).(cellname);
                datasum.(experiment)= calculateelfiz_new(iv.(experiment),data.(experiment));
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
        end
    end
    progressbar(i/length(a));
end
clear a temp loadtemp;

[valuenames,fnames]=getvaluenames(datasum,1);
if ~(getallvalues==1)
    %     for i=1:length(valuenamesidontcare)  %If not valuenames are need %
    %         valuenames(ismember(valuenames,valuenamesidontcare(i)))=[];
    %     end
    valuenames=valuenamesineed;
end
[out,datasum]=getdatafromdatasum_sz(fnames, datasum, valuenames);

[nrows,ncols]=size(out);

outfilename='alldata.dat';
outfile=fopen(outfilename, 'w');

for col=1:ncols-1,
    fprintf(outfile, '%s\t', out{1,col});
end
fprintf(outfile, '%s\n', out{1,ncols});

for row=2:nrows,
    fprintf(outfile, '%s\t', out{row,1});
    for col=2:ncols-1,
        fprintf(outfile, '%f\t', out{row,col});
    end
    fprintf(outfile, '%f\n', out{row,ncols});
end

fclose(outfile);
