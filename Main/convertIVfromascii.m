%clear all
close all
%%%%%%%%%%%%
% Don't forget to add folders to path!!!
% Important to add a subfolder in IV folder
%%%%%%%%%%%%
%addpath('C:\Users\Oliver\Downloads\mobile_taxonomy_170410');
addpath('C:\Users\Oliver\Desktop\Egyetem\7\Szakdoga');
%asciipath='/media/gtamas/Elements/HEKAdata/gulyas';
%asciipath=uigetdir('','Choose directory of ASCII files');
%savepath='/media/gtamas/Elements/MATLABdata/IV/gulyas';
%savepath=uigetdir('','Choose directory of MATLAB files (MATLABdata/IV/*example*)');
%sourcepath=uigetdir('','Choose MATLABdata directory (eg. MATLABdata_allhippo)');
% protocols = GetGoogleSpreadsheet('0Ai5nOoq3uFrCdDFHMElrX0ZRQkNpQWhkbWZQbFk2OUE');
% expnames=protocols(:,find(strcmp(protocols(1,:),'Experiment')));
% stepnames=protocols(:,find(strcmp(protocols(1,:),'protocol')));
%asciipath = 'C:\Users\Oliver\Downloads\mobile_taxonomy_170410\mobile_taxonomy\HEKAdata\zsolt\PV_step_default\Selected';
%savepath='C:\Users\Oliver\Downloads\mobile_taxonomy_170410\mobile_taxonomy\MATLABdata_Zsolt_Newest\IV\PV_step_default';
%sourcepath='C:\Users\Oliver\Downloads\mobile_taxonomy_170410\mobile_taxonomy\MATLABdata_Zsolt_Newest';


% asciipath = 'C:\Users\Oliver\Desktop\Egyetem\7\Szakdoga\HBP_data_all_grouped_szakdoga\PV\STEP HBP';
% savepath='C:\Users\Oliver\Desktop\Egyetem\7\Szakdoga\tmp\IV\OTHER\STEP_HBP';
% sourcepath='C:\Users\Oliver\Desktop\Egyetem\7\Szakdoga\tmp';
%asciipath = 'C:\Users\Oliver\Desktop\Egyetem\7\Szakdoga\Elphys_all_fixed\Other\STEP_LONG';
%savepath='C:\Users\Oliver\Desktop\Egyetem\7\Szakdoga\Elphys_all_fixed_data_again\IV\Other'; %Important to add a subfolder in IV folder
%sourcepath='C:\Users\Oliver\Desktop\Egyetem\7\Szakdoga\Elphys_all_fixed_data_again';

asciipath = 'C:\Users\Oliver\Desktop\Egyetem\7\Szakdoga\Makara_Judit_rat\Makara_Judit_rat\CA1';
savepath='C:\Users\Oliver\Desktop\Egyetem\7\Szakdoga\Makara\IV\CA1'; %Important to add a subfolder in IV folder
sourcepath='C:\Users\Oliver\Desktop\Egyetem\7\Szakdoga\Makara';
% STEP LONG:
%current=[+10,-10,20,-20,30,-30,40,-40,50,-50,60,-60,70,-70,80,-80,90,-90,100,-100,150,200,250,300,400,500,600];
%STEP HBP:
%current=[+20,-20,40,-40,60,-60,80,-80,100,-100,120,-120,140,-140,160,-160,180,-180,200,220,240,260,280,300,400,500,600];
%
current = [0.05, 0.1, 0.15, 0.2, 0.25, -0.05, -0.1, -0.15, -0.2, -0.25];
sampleinterval = 0.0002;	

Folder_Checker_Taxonomy(sourcepath,savepath);

temp=dir(asciipath);
fnames=[];
for i=1:length(temp)
    if ~temp(i).isdir
        fnames{length(fnames)+1}=temp(i).name;
    end
end
lastfname='none';

for i=1:length(fnames)
    fname=char(fnames(i));
    fnamefull=char(fnames(i));
%     if any(strcmp(expnames,fnamefull)) & ~isempty(strfind(char(stepnames(find(strcmp(expnames,fnamefull)))),'Small'))
       
%     else
%         current=stepcurrent;
%     end
    temp=strfind(fname,'.txt');
    if isempty(temp)
        temp=strfind(fname,'_columnwise.txt');
    end
    
    fname=fname(1:temp-1);
    
    
    temp=strfind(fname,'c');
    cellname=['cell_',fname(temp+1:end)];
    
    fname=fname(1:temp-1);
    if exist('iv','var') && strcmp(fname,lastfname)
    elseif exist('iv','var')
        save([savepath,filesep,lastfname],'iv');
        clear iv
    end
    lastfname=fname;
    clear temp;
    temp=dlmread([asciipath,filesep,fnamefull]);
    temp=temp/1000;
%     while mean(temp(1,:))<-0.3
%         temp=temp/10;
%     end
    current=current(1:size(temp,2));
    
    
    iv.(cellname).current = current;
    iv.(cellname).timertime=iv.(cellname).current*0;
    currstodelete=[];
    currdeleteimmediately=[];
    
    % Rough preselection of invalid data columns (current step voltages)
%     for ii=1:size(temp,2)
%         if min(temp(:,currentorder(ii))) < -0.12
%             currstodelete(length(currstodelete)+1)=ii;
%         elseif sum(temp(:,currentorder(ii))==0)>1000
%             currdeleteimmediately(length(currdeleteimmediately)+1)=ii;
%         end
%     end
    
    if ~isempty(currdeleteimmediately)
        currentorder(currdeleteimmediately)=[];
    end
    if ~isempty(currstodelete)
        currentorder(1:currstodelete(end))=[];
    end
    for ii=1:length(current)
        iv.(cellname).(['v',num2str(ii)])=temp(:,ii);
    end
    
    if ~isempty(currstodelete)
        iv.(cellname).current=current(current);
    end
    iv.(cellname).realcurrent=iv.(cellname).current;
    iv.(cellname).bridged=1;
    iv.(cellname).holding=0;
    iv.(cellname).sweepnum=ii;
    iv.(cellname).time=[0:sampleinterval:sampleinterval*(size(temp,1)-1)];
    iv.(cellname).time=iv.(cellname).time';
    iv.(cellname).segment=[100,300,400];
end
spath = [savepath filesep lastfname];
save(spath, 'iv');
clear iv

%sourcepath='/home/halasz/mobile_taxonomy/MATLABdata_Katona';
%sourcepath='/home/heivi/mobile_taxonomy/MATLABdata_zsolt_new';
HEKA_getgoodIVs_main(sourcepath);
% 
