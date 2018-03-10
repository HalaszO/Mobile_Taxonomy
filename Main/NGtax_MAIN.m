clear all
close all
%%%%%%%%%%%%
% Don't forget to add folders to path!
%%%%%%%%%%%%


dothedatasum = 0; %conditional variables
dotheanalysis = 1; 
%addpath('C:\Users\Oliver\Desktop\Egyetem\7\Szakdoga');
% Setting source paths
sourcepath = uigetdir('','Location of MATLABdata library');
ivpath = [sourcepath, filesep 'IV'];
datapath = [sourcepath,filesep 'data'];
picturepath_iv = [sourcepath, filesep 'pictures' filesep 'IV'];
%logfile = [sourcepath, filesep 'Log.meta'];


% Setting some basic values
minsweepnum = 6; % Min number of sweeps required
v0max = -.02; % Deliberately chosen possible max value of the baseline
holdingmin = -500; 
holdingmax = 5000;
disp(sourcepath)
metadata_dir = [sourcepath, filesep 'metadata'];
if ~exist(metadata_dir,'dir') 
    mkdir(metadata_dir); 
end
cd(metadata_dir);
disp(sourcepath)
temp = dir;
metafiles = [];
%logFile = fopen(logfile,'w');

% Checks if contents of the directory are files or subfolders
% If file then append to metafiles
for tempi = 1:length(temp)
    if ~temp(tempi).isdir;
        metafiles{length(metafiles) + 1} = temp(tempi).name;
    end
end
disp(metafiles)

% Going through metafiles 
for metafilenum = 1:length(metafiles)
    % Calling progressbar
    progressbar(['going through   ',num2str(length(metafiles)),'  metafiles'],char(metafiles(metafilenum)), ...
    'passive membrane properties','active membrane properties','active membrane properties no.2');
    progressbar((metafilenum-1)/length(metafiles),[],[],[],[]);
    
    GoodIV_struct = load([sourcepath, filesep ,'metadata', filesep ,char(metafiles(metafilenum))]);
    GoodIV_struct = GoodIV_struct.GoodIVs; % Loading GoodIVs from metafiles
    setupnames = fieldnames(GoodIV_struct); 
    for setupnum = 1:length(setupnames) % Going through struct fields within GoodIVs
        setupname = char(setupnames(setupnum));
        i = 0;
        while i < size(GoodIV_struct.(setupname),2) 
            % Going through the file records
            % Record.field contains the values!
            
            i = i+1;            
            if ~isempty(GoodIV_struct.(setupname)(1,i).dir) 
                %&& ~any(strfind(GoodIV_struct.(setupname)(1,i).fname,'ásolat'))
                fpath = GoodIV_struct.(setupname)(1,i).dir; % IV file location
                pathvar = fpath(length(ivpath)+1:end); % e.g. CCK_Population
                fname = GoodIV_struct.(setupname)(1,i).fname; % e.g. s130201.mat
                cellnames = GoodIV_struct.(setupname)(1,i).ivnames; % Cell name in measurement ID (e.g. cell_0101_ch2)
                recdate = fname(1:strfind(fname,'.')-1); % e.g. s130201
                disp(recdate);
                recdate = ['elfiz',recdate];
                recdate(recdate == '-') = '_'; % Change '-' for '_' -> format specific
                temp = dir([datapath,filesep,pathvar,filesep,'data_iv_',fname]);  % matlab metadata regarding the file (bytes, isdir...)
                % 'data_iv_': data file containing ALL raw IVs and features
                
                %if isempty(temp) && dotheanalysis == 1;
                if dotheanalysis == 1; % if data file doesn't exist YET?
                    temp = load([fpath,filesep,fname]); % load IV file from record.dir 
                    iv.(recdate) = temp.iv; % ...
                    clear temp;
                    
%                     tempcom = char(textscan(commfile,'%s'));
                    %if ~strcmp(GoodIV_struct.(setupname)(1,i).fname, tempcom)  % metafile, shows the current file
%                         outfile = fopen(commfile, 'wt' ); 
%                         fprintf(outfile,'%s\r\n', fname);
%                         fclose(outfile);
                        
                        for ii=1:length(cellnames)
                            cellname=char(cellnames(ii));
                            if 0
%                                 isnan(GoodIV_struct.(setupname)(1,i).holding(ii)) || isempty(GoodIV_struct.(setupname)(1,i).v0) || GoodIV_struct.(setupname)(1,i).v0(ii)>v0max ...
%                                     || GoodIV_struct.(setupname)(1,i).holding(ii)<holdingmin || GoodIV_struct.(setupname)(1,i).holding(ii)>holdingmax || ...
%                                     GoodIV_struct.(setupname)(1,i).sweepnum(ii)<minsweepnum || GoodIV_struct.(setupname)(1,i).lastcurrent(ii)<=0
%                                 disp([cellname,' -- fail']);
                            else
                                templength=[];
                                for tempi=1:iv.(recdate).(cellname).sweepnum
                                    templength(tempi)=length(iv.(recdate).(cellname).(['v',num2str(iv.(recdate).(cellname).sweepnum)]));
                                end
                                if ~sum(templength==length(iv.(recdate).(cellname).time))
                                    disp([cellname,' -- fail not equal sweeps']);
                                elseif iv.(recdate).(cellname).time(end)<sum(iv.(recdate).(cellname).segment(1:2)/1000)
                                    disp([cellname,' -- fail not correct segments']);
                                else
                                    
                                    disp(cellname);
                                    
                                    % Important part, calling the functions
                                    
                                    data.(recdate).(cellname).pass=mpassive_sag(iv.(recdate).(cellname),cellname);
                                    data.(recdate).(cellname).HH = mHH(iv.(recdate).(cellname),data.(recdate).(cellname).pass.dvrs,data.(recdate).(cellname).pass.vrs);
                                    %%data.(recdate).(cellname)=offsetvoltage(data.(recdate).(cellname),iv.(recdate).(cellname)); 
                                    %%Incorporated into .mHH
                                    datasum.(recdate).(cellname)= calculateelfiz_new(iv.(recdate).(cellname),data.(recdate).(cellname));
%                                     plotiv(cellname, iv.(recdate), datasum.(recdate), data.(recdate), 1, 1000);
%                                     if ~exist([picturepath_iv,pathvar],'dir') 
%                                         mkdir([picturepath_iv,pathvar]); 
%                                     end
%                                     saveas(gcf,[picturepath_iv,pathvar,filesep,recdate,'_',cellname,'_halol.fig']);
%                                     saveas(gcf,[picturepath_iv,pathvar,filesep,recdate,'_',cellname,'_halol.png']);
%                                    spikestoDelete = data.(recdate).(cellname).HH.NumofDeletedSpikes;
%                                    for sptd = 1:length(spikestoDelete)
 %                                       if ~isempty(spikestoDelete(sptd)) && spikestoDelete(sptd) ~= 0
 %                                           fprintf(logFile, 'From %s%s %d sweep %d spikes were deleted\n',recdate,cellname,sptd,spikestoDelete(sptd));
  %                                      end
  %                                  end
                                end
                            end
                        end
                        if exist('data','var');
                            if ~exist([datapath,pathvar],'dir') 
                                mkdir([datapath,pathvar]);
                            end
                            save([datapath,pathvar,filesep,'data_iv_',fname], 'data', 'iv');
                            clear data
                        end
                        clear iv
                    %end
                end
            end
            %%% Creation of datasum
            if ~exist('datasum','var')
                datasum = struct;
            end
            if dothedatasum == 1 && ~isempty(dir([datapath,pathvar,filesep,'data_iv_',fname])) 
                temp=load([datapath,pathvar,filesep,'data_iv_',fname]);
                if  ~isfield(temp.data,recdate) 
                    delete([datapath,pathvar,filesep,'data_iv_',fname]);
                    i=i-1;
                else
                    cellnames=fieldnames(temp.data.(recdate));
                    for ii=1:length(cellnames)
                        cellname=char(cellnames(ii));
                        if sum(temp.data.(recdate).(cellname).HH.apnum)>0
                            datasum.(recdate).(cellname)= calculateelfiz_new(temp.iv.(recdate).(cellname),temp.data.(recdate).(cellname));
                            datasum.(recdate).(cellname).pathandfname={[datapath,pathvar],[filesep,'data_iv_',fname]};
                            datasum.(recdate).(cellname).timesincemidnight=GoodIV_struct.(setupname)(1,i).timestart(ii);
                            %plotiv(cellname, temp.iv.(recdate), datasum.(recdate), temp.data.(recdate), 1, 1000);
                            %saveas(gcf,[picturepath_iv,pathvar,filesep,recdate,'_',cellname,'_halol.fig']);
                            %saveas(gcf,[picturepath_iv,pathvar,filesep,recdate,'_',cellname,'_halol.jpg']);
                        end
                    end
                end
                clear temp;
            end
            %%% datasum legyártása még1x
            progressbar([],i/size(GoodIV_struct.(setupname),2));
        end
    end
    progressbar(metafilenum/length(metafiles));
end
if dothedatasum==1
    save([datapath,filesep,'datasum_PV'],'datasum');
end

%fclose(logFile);
