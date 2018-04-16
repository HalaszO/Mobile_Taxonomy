function HEKA_getgoodIVs_main(sourcepath)
%% nyers IVk becuppantasa
% sourcepath='/media/gtamas/Elements/MATLABdata';
ivpath=[sourcepath,filesep 'IV'];
analyseallfiles=1;

if ~analyseallfiles==1
    cd(ivpath);
    rawivdirs=uipickfiles;
else
    rawivdirs={ivpath};
end
prevdirnum=length(rawivdirs)-1;
while length(rawivdirs)>prevdirnum
    prevdirnum=length(rawivdirs);
    for i=1:length(rawivdirs) %% megnézzük, hogy hány mappa van a kiválasztott mappákon belül
        rawivdir=char(rawivdirs(i));
            cd(rawivdir);
            temp=dir;
            for j=1:length(temp)
                if temp(j).isdir==1 && ~or(strcmp(temp(j).name,'.'),strcmp(temp(j).name,'..')) && ~any(strcmp([rawivdir, filesep ,temp(j).name],rawivdirs))
                    rawivdirs{length(rawivdirs)+1}=[rawivdir,filesep,temp(j).name];
                end
            end
    end
end
rawivfnames=struct;
for i=1:length(rawivdirs)
    rawivdir=char(rawivdirs(i));
    if ~strcmp(rawivdir,ivpath)
        temp=[rawivdir(length(ivpath)+2:end),filesep];
        setupname=['s_',temp(1:find(temp== filesep,1,'first')-1)];
        if any(strfind(setupname,'-'))
            setupname(setupname=='-')='_';
        end
        if ~isfield(rawivfnames,setupname)
            rawivfnames.(setupname)=[];
        end
        cd(rawivdir);
        temp=dir;
        for j=1:length(temp)
            if temp(j).isdir==0 && strcmp(temp(j).name(end-2:end),'mat')
                tempdb=size(rawivfnames.(setupname),1)+1;
                rawivfnames.(setupname){tempdb,1}=rawivdir;
                rawivfnames.(setupname){tempdb,2}=temp(j).name;
                rawivfnames.(setupname){tempdb,3}=datenum(char(temp(j).date),0);
                if length(temp(j).name)>9
                    rawivfnames.(setupname){tempdb,4}=datenum(temp(j).name(1:6),'yymmdd');
                    rawivfnames.(setupname){tempdb,5}=abs(cell2mat(rawivfnames.(setupname)(tempdb,3))-cell2mat(rawivfnames.(setupname)(tempdb,4)));
                else
                    rawivfnames.(setupname){tempdb,4}='NaN';
                    rawivfnames.(setupname){tempdb,5}=0;
                end
            end
        end
    end
end
clear temp tempdb prevdirnum i j rawivdir rawivdirs analyseallfiles
%% holding, step kiszedése
GoodIVs=struct;
setupnames=fieldnames(rawivfnames);
for setupnum=1:length(setupnames)
    setupname=char(setupnames(setupnum));
    progressbar([setupname,'  :  ',num2str(setupnum/length(setupnames)*100),'%']);
    if ~isfield(GoodIVs,setupname)
        GoodIVs.(setupname)=[];
    end
    for i=1:size(rawivfnames.(setupname),1)
        rawivdir=char(rawivfnames.(setupname)(i,1));
        rawivfname=char(rawivfnames.(setupname)(i,2));
        tempp=dir([rawivdir,filesep,rawivfname]);
        temp=load([rawivdir,filesep,rawivfname]);
        ivnames=fieldnames(temp.iv);
        neededivnames=ivnames;
        ivname=char(ivnames(1));
        %%%% megadjuk a legnagyobb std-vel rendelkező csatornák nevét
%         while ~isempty(ivnames)
%             
%             group=str2double(ivname(strfind(ivname,'g')+1:strfind(ivname,'_s')-1));
%             series=str2double(ivname(strfind(ivname,'_s')+2:strfind(ivname,'_c')-1));
%             channel=str2double(ivname(strfind(ivname,'_')+2:end));
%             sameseries=ivnames(strncmp(ivname,ivnames,end));
%             vstd=[];
%             for tempi=1:length(sameseries)
%                 ivname=char(sameseries(tempi));
%                 tempvstd=[];
%                 for tempj=1:temp.iv.(ivname).sweepnum;
%                     if isfield(temp.iv.(ivname),(['v',num2str(tempj)]))
%                         tempvstd(tempj)=std(temp.iv.(ivname).(['v',num2str(tempj)]));
%                     end
%                 end
%                 vstd(tempi)=mean(tempvstd);
%                 ivnames(1)=[];
%             end
%             if ~isnan(max(vstd))
%                 neededivnames{length(neededivnames)+1}=char(sameseries(find(vstd==max(vstd),1,'first')));
%             end
%             
%             
%         end
        %%%% megadjuk a legnagyobb std-vel rendelkező csatornák nevét
        if ~isempty(neededivnames)
            GoodIVs.(setupname)(i).ivfile.lastmodify=tempp.date;
            GoodIVs.(setupname)(i).ivfile.bytes=tempp.bytes;
            GoodIVs.(setupname)(i).fname=rawivfname;
            GoodIVs.(setupname)(i).dir=rawivdir;
            GoodIVs.(setupname)(i).ivnames=neededivnames;
            
            GoodIVs.(setupname)(i).firstcurrent=zeros(length(neededivnames),1);
            GoodIVs.(setupname)(i).lastcurrent=zeros(length(neededivnames),1);
            GoodIVs.(setupname)(i).holding=zeros(length(neededivnames),1);
            GoodIVs.(setupname)(i).step=zeros(length(neededivnames),1);
            GoodIVs.(setupname)(i).sweepnum=zeros(length(neededivnames),1);
            for tempi=1:length(neededivnames)
                neededivname=char(neededivnames(tempi));
                GoodIVs.(setupname)(i).firstcurrent(tempi)=temp.iv.(neededivname).current(1);
                GoodIVs.(setupname)(i).lastcurrent(tempi)=temp.iv.(neededivname).current(end);
                GoodIVs.(setupname)(i).holding(tempi)=temp.iv.(neededivname).holding;
                GoodIVs.(setupname)(i).step(tempi)=temp.iv.(neededivname).current(2)-temp.iv.(neededivname).current(1);
                GoodIVs.(setupname)(i).sweepnum(tempi)=temp.iv.(neededivname).sweepnum;
                if isfield(temp.iv.(ivname),'v1')
                    GoodIVs.(setupname)(i).v0(tempi)=mean(temp.iv.(neededivname).v1(1:find(temp.iv.(neededivname).time>temp.iv.(neededivname).segment(1)/1000,1,'first')));
                end
                GoodIVs.(setupname)(i).sampleinterval(tempi)=temp.iv.(neededivname).time(2)-temp.iv.(neededivname).time(1);
                GoodIVs.(setupname)(i).timestart(tempi)=temp.iv.(neededivname).timertime(1);
                GoodIVs.(setupname)(i).timeend(tempi)=temp.iv.(neededivname).timertime(end);
            end
        end
        progressbar(i/size(rawivfnames.(setupname),1));
    end
        
    save([sourcepath,filesep 'metadata' filesep ,setupname,'.mat'],'GoodIVs');
    GoodIVs=rmfield(GoodIVs,setupname);
end

%%

% %% ellenőrzés
% cd([ivpath,'/metadata/']);
% temp=dir;
% minsweepnum=5;
% v0max=-.03;
% holdingmin=-500;
% holdingmax=5000;
% plotbadones=0;
% badfiles=[];
% badnum=0;
% for i=1:size(GoodIVs,2)
%     badones=0;
%     for j=1:length(GoodIVs(1,i).sweepnum)
%         ivname=char(GoodIVs(1,i).ivnames(j));
%         temp=[];
%         if isempty(GoodIVs(1,i).v0) | GoodIVs(1,i).v0(j)>v0max | GoodIVs(1,i).holding(j)<holdingmin | GoodIVs(1,i).holding(j)>holdingmax | GoodIVs(1,i).sweepnum(j)<minsweepnum
%             badones=badones+1;
%             if isempty(temp)
%                 temp=load([GoodIVs(1,i).dir,'/',GoodIVs(1,i).fname]);
%             end
%             if plotbadones==1
%                 figure(1);
%                 cla;
%                 hold on;
%                 plot(temp.iv.(ivname).time,temp.iv.(ivname).v1);
%                 plot(temp.iv.(ivname).time,temp.iv.(ivname).(['v',num2str(GoodIVs(1,i).sweepnum(j))]));
%                 %[GoodIVs(1,i).v0(j), GoodIVs(1,i).holding(j), GoodIVs(1,i).sweepnum(j)]
%                 pause(3);
%             end
%             badnum=badnum+1;
%         end
%     end
%     if badones==j
%         nextrow=size(badfiles,1)+1;
%         badfiles{nextrow,1}=GoodIVs(1,i).dir;
%         badfiles{nextrow,2}=GoodIVs(1,i).fname;
%     end
%     progressbar(i/size(GoodIVs,2));
% end
% % % % % % % % % %% delete bad files
% % % % % % % % % for i=1:size(badfiles,1)
% % % % % % % % %     [char(badfiles(i,1)),'/',char(badfiles(i,2))]
% % % % % % % % %     delete([char(badfiles(i,1)),'/',char(badfiles(i,2))]);
% % % % % % % % % end