function [datasum] = calculateelfiz_new(iv,data)
%neededsteadyapnum=8;
RS=0;
RScount=0;
datasum.burstspikes=zeros(1,iv.sweepnum);
for i=1:iv.sweepnum
    if -40< iv.current(i) && iv.current(i)<40 % Averaging of RS
    else
        RS=RS+data.pass.rs(i);
        RScount=RScount+1;
    end
    if  data.HH.apnum(i)>1 & isfield(data.HH.(['sweep',num2str(i)]),'interspike') & size(data.HH.(['sweep',num2str(i)]).interspike)>=2 & data.HH.(['sweep',num2str(i)]).apmaxtime(1)<iv.segment(1)/1000+.15
        % Reminder: the 1st element of interspike is the mean, hence the
        % array length matches that of the AP arrays
        if data.HH.(['sweep',num2str(i)]).interspike(2)<70
            datasum.burstspikes(i)=2;
            while datasum.burstspikes(i)+1<=data.HH.apnum(i) && or(data.HH.(['sweep',num2str(i)]).interspike(datasum.burstspikes(i)+1)<data.HH.(['sweep',num2str(i)]).interspike(datasum.burstspikes(i))*2, data.HH.(['sweep',num2str(i)]).interspike(datasum.burstspikes(i)+1)<20)
                datasum.burstspikes(i)=datasum.burstspikes(i)+1;
            end
            restisimin=100;
            if data.HH.apnum(i)>datasum.burstspikes(i)+1
                restisimin=min(data.HH.(['sweep',num2str(i)]).interspike(datasum.burstspikes(i)+2:end));
            end
            if datasum.burstspikes(i)>8 || data.HH.(['sweep',num2str(i)]).apmaxtime(datasum.burstspikes(i))-data.HH.(['sweep',num2str(i)]).apmaxtime(1)>.08 || restisimin<max(data.HH.(['sweep',num2str(i)]).interspike(2:datasum.burstspikes(i)));
                datasum.burstspikes(i)=0;
            end
        end
    end
    
end
datasum.noiselevel=data.pass.noiselevel;
datasum.filterednoiselevel=data.pass.filterednoiselevel;
datasum.sampleinterval=1/data.pass.samplingrate*1000000;
datasum.RS=RS/RScount;
datasum.sagsold=data.pass.rsag./data.pass.rin_old(1:length(data.pass.rsag));
datasum.sags = data.pass.rsag/data.pass.rin_new;
datasum.firstsag=mean(datasum.sags(1:1));
for i=2:4
    if length(datasum.sags)>=i && std(datasum.sags(1:i))<.1
        datasum.firstsag=mean(datasum.sags(1:i));
    end
end
datasum.rebounds=data.pass.dvrebound./data.pass.dvin(1:size(data.pass.dvrebound,2));
datasum.rebound=datasum.rebounds(1);
if isfield(data.HH,'rhump') %%mert ilyen is el�fordul, ha nincs pozit�v �raml�pcs�, �s r�gt�n t�zel
    datasum.humps=data.HH.rhump./data.pass.rin_old(1:length(data.HH.rhump));
    datasum.HUMP=(datasum.humps(end));
else
    datasum.humps=0;
    datasum.HUMP=0;
end
datasum.Rinolds=data.pass.rin_old(1:5);
datasum.Rin_new=mean(data.pass.rin_new);
datasum.Rinstd=std(datasum.Rinolds(1:3));
% if datasum.Rinstd>datasum.Rin/10 % In case of a large std deviance, we only use the first one
%     datasum.Rin=mean(datasum.Rins(1:1));
%     datasum.Rinstd=0;
% end

datasum.vresting=data.pass.v0*1000;
datasum.v0=mean(data.pass.v0s);
datasum.v0std=std(data.pass.v0s);
datasum.apnum=data.HH.apnum;
[datasum.maxapnum, datasum.maxapsweep]=max(datasum.apnum);
datasum.maxapcurrent=iv.realcurrent(datasum.maxapsweep);
datasum.reobasesweep=data.HH.reobasesweep;
if datasum.reobasesweep>0
    datasum.reobase=iv.realcurrent(datasum.reobasesweep);%datasum.reobasesweep*20-120;
    lastzerosweep=max(datasum.reobasesweep-1,1);
    datasum.fIslope=datasum.maxapnum/(datasum.maxapcurrent-iv.realcurrent(lastzerosweep));
    datasum.apreduce=1-datasum.apnum(end)/datasum.maxapnum;
else
    datasum.reobase=NaN;
end


if length(data.pass.tau_new)>=5
datasum.tau_news=data.pass.tau_new(1:5);
else
    datasum.tau_news=data.pass.tau_new;
end
datasum.tauold1=data.pass.tauold(1);
datasum.taunoldfail1=data.pass.tauoldfail(1);

if any(data.HH.apnum>0)
    
    datasum.apnumfromreobase=datasum.apnum(datasum.reobasesweep:end);
    datasum.steadysweep=data.HH.steadysweep;
    datasum.steadyAPnum=datasum.apnum(datasum.steadysweep);
    datasum.goodsteadyAPnum=datasum.steadyAPnum-datasum.burstspikes(datasum.steadysweep)-1;
    
    
    %%%%%%%%%%%%%%%%%%AP UT�NI ESEM�NYEK!!!!!!!!!!!!!!
    %%%%%%%%%%%%%%%%%%AP UT�NI ESEM�NYEK!!!!!!!!!!!!!!
    mindiff=0.0005; %ami alatt nem vesz�nk m�r k�l�nbs�get az AP ut�ni esem�nyekre
    %minVdiff=.2/1000; %amilyen fesz�lts�gk�l�nbs�g alatt m�r nem vesz�nk k�l�nbs�get..
    if ~isempty(datasum.reobasesweep)
        datasum.threshold=data.HH.(['sweep',num2str(datasum.reobasesweep)]).thresh_corrected(1)*1000;
        datasum.apampl=data.HH.(['sweep',num2str(datasum.reobasesweep)]).apamplitude(1);
        datasum.aprisetime=data.HH.(['sweep',num2str(datasum.reobasesweep)]).apmaxtime(1)-data.HH.(['sweep',num2str(datasum.reobasesweep)]).threshtime(1);
        datasum.aphw=data.HH.(['sweep',num2str(datasum.reobasesweep)]).halfwidth(1);
        
        datasum.apwidth=(data.HH.(['sweep',num2str(datasum.reobasesweep)]).apendtime(1)-data.HH.(['sweep',num2str(datasum.reobasesweep)]).threshtime(1))*1000;
        datasum.apstartenddiff=(data.HH.(['sweep',num2str(datasum.reobasesweep)]).thresh(1)-data.HH.(['sweep',num2str(datasum.reobasesweep)]).apend(1))*1000;
        datasum.dvmax=data.HH.(['sweep',num2str(datasum.reobasesweep)]).dvmax(1);
        datasum.dvmaxv=data.HH.(['sweep',num2str(datasum.reobasesweep)]).dvmaxv_corrected(1)*1000;
        datasum.dvmin=data.HH.(['sweep',num2str(datasum.reobasesweep)]).dvmin(1);
        datasum.dvminv=data.HH.(['sweep',num2str(datasum.reobasesweep)]).dvminv_corrected(1)*1000;
        
        if (data.HH.(['sweep',num2str(datasum.reobasesweep)]).apmaxtime(1)-data.HH.(['sweep',num2str(datasum.reobasesweep)]).ahptime(1))>mindiff %% We only care about AHP-ADP if the next threshold is not reached
            datasum.ahpampl=(data.HH.(['sweep',num2str(datasum.reobasesweep)]).apend(1)-data.HH.(['sweep',num2str(datasum.reobasesweep)]).ahpv(1))*1000;
            datasum.ahpwidth=(data.HH.(['sweep',num2str(datasum.reobasesweep)]).ahptime(1)-data.HH.(['sweep',num2str(datasum.reobasesweep)]).apendtime(1))*1000;
            if (data.HH.(['sweep',num2str(datasum.reobasesweep)]).maxtime(1)-data.HH.(['sweep',num2str(datasum.reobasesweep)]).adptime(1))>mindiff
                datasum.adpampl=(data.HH.(['sweep',num2str(datasum.reobasesweep)]).adpv(1)-data.HH.(['sweep',num2str(datasum.reobasesweep)]).ahpv(1))*1000;
                datasum.adpwidth=(data.HH.(['sweep',num2str(datasum.reobasesweep)]).adptime(1)-data.HH.(['sweep',num2str(datasum.reobasesweep)]).ahptime(1))*1000;
                if exist(strcat('data.HH.([sweep',num2str(datasum.reobasesweep),']).ahptimeslow'), 'var') && (data.HH.(['sweep',num2str(datasum.reobasesweep)]).maxtime(1)-data.HH.(['sweep',num2str(datasum.reobasesweep)]).ahptimeslow(1))>mindiff
                    datasum.ahpslowampl=(data.HH.(['sweep',num2str(datasum.reobasesweep)]).adpv(1)-data.HH.(['sweep',num2str(datasum.reobasesweep)]).ahpvslow(1))*1000;
                    datasum.ahpslowwidth=(data.HH.(['sweep',num2str(datasum.reobasesweep)]).ahptimeslow(1)-data.HH.(['sweep',num2str(datasum.reobasesweep)]).adptime(1))*1000;
                else
                    datasum.ahpslowampl=0;
                    datasum.ahpslowwidth=0;
                end
            else
                datasum.adpampl=0;
                datasum.adpwidth=0;
                datasum.ahpslowampl=0;
                datasum.ahpslowwidth=0;
            end
        else
            datasum.ahpampl=0;
            datasum.ahpwidth=0;
            datasum.adpampl=0;
            datasum.adpwidth=0;
            datasum.ahpslowampl=0;
            datasum.ahpslowwidth=0;
        end
        if (data.HH.(['sweep',num2str(datasum.reobasesweep)]).adptime(1)-data.HH.(['sweep',num2str(datasum.reobasesweep)]).ahptime(1))<mindiff
            datasum.adpampl=0;
            datasum.adpwidth=0;
        end
        if exist(strcat('data.HH.([sweep',num2str(datasum.reobasesweep),']).ahptimeslow'), 'var') && (data.HH.(['sweep',num2str(datasum.reobasesweep)]).ahptimeslow(1)-data.HH.(['sweep',num2str(datasum.reobasesweep)]).adptime(1))<mindiff
            datasum.ahpslowampl=0;
            datasum.ahpslowwidth=0;
        end
        if datasum.adpwidth>20 %%%%%%%%%% inplausable width
            datasum.adpampl=0;
            datasum.adpwidth=0;
            if exist(strcat('data.HH.([sweep',num2str(datasum.reobasesweep),']).ahptimeslow'), 'var')
                datasum.ahpslowampl=0;
                datasum.ahpslowwidth=0;
            end
        end
        datasum.firstapmaxtime=data.HH.(['sweep',num2str(datasum.reobasesweep)]).apmaxtime(1)*1000-iv.segment(1);
        datasum.burstinterval=NaN;
        if datasum.burstspikes(datasum.steadysweep)>0
            datasum.burstinterval=mean(data.HH.(['sweep',num2str(datasum.steadysweep)]).interspike(2:datasum.burstspikes(datasum.steadysweep)));
        end
        datasum.steadyinterval=100;
        if datasum.apnum(datasum.steadysweep)>datasum.burstspikes(datasum.steadysweep)+1
            datasum.steadyinterval=mean(data.HH.(['sweep',num2str(datasum.steadysweep)]).interspike(datasum.burstspikes(datasum.steadysweep)+2:end));
        end
        datasum.compfail=mean(data.HH.(['sweep',num2str(datasum.steadysweep)]).compfail);
        if datasum.compfail>0
            datasum.compfail=1;
        end
    end
    
    
    
    
    if ~isempty(datasum.steadysweep) %%%%%%%%%%%%%%%%%%%%%%steady state APs
        for i=1:datasum.steadyAPnum
            
            datasum.steadyaphws(i)=data.HH.(['sweep',num2str(datasum.steadysweep)]).halfwidth(i);
            datasum.steadyapampls(i)=data.HH.(['sweep',num2str(datasum.steadysweep)]).apamplitude(i);
            datasum.steadythresholds(i)=data.HH.(['sweep',num2str(datasum.steadysweep)]).thresh_corrected(i)*1000;
            datasum.steadyapwidths(i)=(data.HH.(['sweep',num2str(datasum.steadysweep)]).apendtime(i)-data.HH.(['sweep',num2str(datasum.steadysweep)]).threshtime(i))*1000;
            datasum.steadyapstartenddiffs(i)=(data.HH.(['sweep',num2str(datasum.steadysweep)]).thresh(i)-data.HH.(['sweep',num2str(datasum.steadysweep)]).apend(i))*1000;
            
            if (data.HH.(['sweep',num2str(datasum.steadysweep)]).maxtime(i)-data.HH.(['sweep',num2str(datasum.steadysweep)]).ahptime(i))>mindiff %% Csak akkor foglalkozunk AHP-ADP-vel, ha nem �rte el a k�vetkez� thresholdot
                datasum.steadyahpampls(i)=(data.HH.(['sweep',num2str(datasum.steadysweep)]).apend(i)-data.HH.(['sweep',num2str(datasum.steadysweep)]).ahpv(i))*1000;
                datasum.steadyahpwidths(i)=(data.HH.(['sweep',num2str(datasum.steadysweep)]).ahptime(i)-data.HH.(['sweep',num2str(datasum.steadysweep)]).apendtime(i))*1000;
                if (data.HH.(['sweep',num2str(datasum.steadysweep)]).maxtime(i)-data.HH.(['sweep',num2str(datasum.steadysweep)]).adptime(i))>mindiff
                    datasum.steadyadpampls(i)=(data.HH.(['sweep',num2str(datasum.steadysweep)]).adpv(i)-data.HH.(['sweep',num2str(datasum.steadysweep)]).ahpv(i))*1000;
                    datasum.steadyadpwidths(i)=(data.HH.(['sweep',num2str(datasum.steadysweep)]).adptime(i)-data.HH.(['sweep',num2str(datasum.steadysweep)]).ahptime(i))*1000;
                    if (data.HH.(['sweep',num2str(datasum.steadysweep)]).maxtime(i)-data.HH.(['sweep',num2str(datasum.steadysweep)]).ahptimeslow(i))>mindiff
                        datasum.steadyahpslowampls(i)=(data.HH.(['sweep',num2str(datasum.steadysweep)]).adpv(i)-data.HH.(['sweep',num2str(datasum.steadysweep)]).ahpvslow(i))*1000;
                        datasum.steadyahpslowwidths(i)=(data.HH.(['sweep',num2str(datasum.steadysweep)]).ahptimeslow(i)-data.HH.(['sweep',num2str(datasum.steadysweep)]).adptime(i))*1000;
                    else
                        datasum.steadyahpslowampls(i)=0;
                        datasum.steadyahpslowwidths(i)=0;
                    end
                else
                    datasum.steadyadpampls(i)=0;
                    datasum.steadyadpwidths(i)=0;
                    datasum.steadyahpslowampls(i)=0;
                    datasum.steadyahpslowwidths(i)=0;
                end
            else
                datasum.steadyahpampls(i)=0;
                datasum.steadyahpwidths(i)=0;
                datasum.steadyadpampls(i)=0;
                datasum.steadyadpwidths(i)=0;
                datasum.steadyahpslowampls(i)=0;
                datasum.steadyahpslowwidths(i)=0;
            end
            
            if abs(data.HH.(['sweep',num2str(datasum.steadysweep)]).adptime(i)-data.HH.(['sweep',num2str(datasum.steadysweep)]).ahptime(i))<mindiff
                datasum.steadyadpampls(i)=0;
                datasum.steadyadpwidths(i)=0;
            end
            if abs(data.HH.(['sweep',num2str(datasum.steadysweep)]).ahptimeslow(i)-data.HH.(['sweep',num2str(datasum.steadysweep)]).adptime(i))<mindiff
                datasum.steadyahpslowampls(i)=0;
                datasum.steadyahpslowwidths(i)=0;
            end
            if datasum.steadyadpwidths(i)>20 %%%%%%%%%% az irre�lisan t�voli ADP eset�n kinull�zok mindent ut�na
                datasum.steadyadpampls(i)=0;
                datasum.steadyadpwidths(i)=0;
                datasum.steadyahpslowampls(i)=0;
                datasum.steadyahpslowwidths(i)=0;
            end
        end
        
        if datasum.burstspikes(datasum.steadysweep)+2<datasum.apnum(datasum.steadysweep)
            datasum.steadyaphw = mean(deleteoutliers(datasum.steadyaphws(datasum.burstspikes(datasum.steadysweep)+1:end-1)));
            datasum.steadyapampl = mean(deleteoutliers(datasum.steadyapampls(datasum.burstspikes(datasum.steadysweep)+1:end-1)));
            datasum.steadythreshold = mean(deleteoutliers(datasum.steadythresholds(datasum.burstspikes(datasum.steadysweep)+1:end-1)));
            datasum.steadythresholddiff = mean(deleteoutliers(datasum.steadythresholds(datasum.burstspikes(datasum.steadysweep)+1:end-1)-data.pass.v0s(datasum.steadysweep)*1000));
            datasum.steadyapwidth = mean(deleteoutliers(datasum.steadyapwidths(datasum.burstspikes(datasum.steadysweep)+1:end-1)));
            datasum.steadyapstartenddiff = mean(deleteoutliers(datasum.steadyapstartenddiffs(datasum.burstspikes(datasum.steadysweep)+1:end-1)));
            datasum.steadyahpampl = mean(deleteoutliers(datasum.steadyahpampls(datasum.burstspikes(datasum.steadysweep)+1:end-1)));
            datasum.steadyahpwidth = mean(deleteoutliers(datasum.steadyahpwidths(datasum.burstspikes(datasum.steadysweep)+1:end-1)));
            datasum.steadyadpampl = mean(deleteoutliers(datasum.steadyadpampls(datasum.burstspikes(datasum.steadysweep)+1:end-1)));
            datasum.steadyadpwidth = mean(deleteoutliers(datasum.steadyadpwidths(datasum.burstspikes(datasum.steadysweep)+1:end-1)));
            datasum.steadyahpslowampl = mean(deleteoutliers(datasum.steadyahpslowampls(datasum.burstspikes(datasum.steadysweep)+1:end-1)));
            datasum.steadyahpslowwidth = mean(deleteoutliers(datasum.steadyahpslowwidths(datasum.burstspikes(datasum.steadysweep)+1:end-1)));
        elseif datasum.burstspikes(datasum.steadysweep)==datasum.apnum(datasum.steadysweep)
            datasum.steadyaphw = datasum.steadyaphws(end);
            datasum.steadyapampl = datasum.steadyapampls(end);
            datasum.steadythreshold = datasum.steadythresholds(end);
            datasum.steadythresholddiff = datasum.steadythresholds(end)-data.pass.v0s(datasum.steadysweep)*1000;
            datasum.steadyapwidth = datasum.steadyapwidths(end);
            datasum.steadyapstartenddiff = datasum.steadyapstartenddiffs(end);
            datasum.steadyahpampl = datasum.steadyahpampls(end);
            datasum.steadyahpwidth = datasum.steadyahpwidths(end);
            datasum.steadyadpampl = datasum.steadyadpampls(end);
            datasum.steadyadpwidth = datasum.steadyadpwidths(end);
            datasum.steadyahpslowampl = datasum.steadyahpslowampls(end);
            datasum.steadyahpslowwidth = datasum.steadyahpslowwidths(end);
        else
            datasum.steadyaphw = datasum.steadyaphws(datasum.burstspikes(datasum.steadysweep)+1);
            datasum.steadyapampl = datasum.steadyapampls(datasum.burstspikes(datasum.steadysweep)+1);
            datasum.steadythreshold = datasum.steadythresholds(datasum.burstspikes(datasum.steadysweep)+1);
            datasum.steadythresholddiff= datasum.steadythresholds(datasum.burstspikes(datasum.steadysweep)+1)-data.pass.v0s(datasum.steadysweep)*1000;
            datasum.steadyapwidth = datasum.steadyapwidths(datasum.burstspikes(datasum.steadysweep)+1);
            datasum.steadyapstartenddiff = datasum.steadyapstartenddiffs(datasum.burstspikes(datasum.steadysweep)+1);
            datasum.steadyahpampl = datasum.steadyahpampls(datasum.burstspikes(datasum.steadysweep)+1);
            datasum.steadyahpwidth = datasum.steadyahpwidths(datasum.burstspikes(datasum.steadysweep)+1);
            datasum.steadyadpampl = datasum.steadyadpampls(datasum.burstspikes(datasum.steadysweep)+1);
            datasum.steadyadpwidth = datasum.steadyadpwidths(datasum.burstspikes(datasum.steadysweep)+1);
            datasum.steadyahpslowampl = datasum.steadyahpslowampls(datasum.burstspikes(datasum.steadysweep)+1);
            datasum.steadyahpslowwidth = datasum.steadyahpslowwidths(datasum.burstspikes(datasum.steadysweep)+1);
        end
        datasum.ahpwidthnew=datasum.apwidth+datasum.ahpwidth;
        datasum.ahpamplnew=datasum.apstartenddiff+datasum.ahpampl;
        %%%%%%%%%%%%%%%%%%%%%%%FIRING PATTERN!!!!!!!!!!!!!
        interspike=data.HH.(['sweep',num2str(datasum.steadysweep)]).interspike;
        datasum.accomodation=0;
        datasum.ISIchange=0;
        datasum.ISIchangevsapnum=0;
        %     datasum.ISIchangeerror=0;
        if datasum.apnum(datasum.steadysweep)-datasum.burstspikes(datasum.steadysweep)>3
            datasum.accomodation=interspike(end)/interspike(datasum.burstspikes(datasum.steadysweep)+2); % az els� nem burst spike ut�ni ISI �s az ut�ls� ISI ar�nya
            [p] = polyfit(data.HH.(['sweep',num2str(datasum.steadysweep)]).apmaxtime(datasum.burstspikes(datasum.steadysweep)+2:end),data.HH.(['sweep',num2str(datasum.steadysweep)]).interspike(datasum.burstspikes(datasum.steadysweep)+2:end),1);
            %         datasum.ISIchangeerror=mean(abs(data.HH.(['sweep',num2str(datasum.steadysweep)]).apmaxtime(2:end).*p(1)+p(2)-data.HH.(['sweep',num2str(datasum.steadysweep)]).interspike(2:end)));
            datasum.ISIchange=p(1);
            [p] = polyfit(1:length(data.HH.(['sweep',num2str(datasum.steadysweep)]).interspike(datasum.burstspikes(datasum.steadysweep)+2:end)),data.HH.(['sweep',num2str(datasum.steadysweep)]).interspike(datasum.burstspikes(datasum.steadysweep)+2:end),1);
            datasum.ISIchangevsapnum=p(1);
            %         figure
            %         plot(data.HH.(['sweep',num2str(datasum.steadysweep)]).apmaxtime(datasum.burstspikes(datasum.steadysweep)+2:end), data.HH.(['sweep',num2str(datasum.steadysweep)]).interspike(datasum.burstspikes(datasum.steadysweep)+2:end),'.')
            %         hold on
            %         plot(data.HH.(['sweep',num2str(datasum.steadysweep)]).apmaxtime(datasum.burstspikes(datasum.steadysweep)+2:end), data.HH.(['sweep',num2str(datasum.steadysweep)]).apmaxtime(datasum.burstspikes(datasum.steadysweep)+2:end).*p(1)+p(2));
            %         plot(data.HH.(['sweep',num2str(datasum.steadysweep)]).apmaxtime(datasum.burstspikes(datasum.steadysweep)+2:end),datasum.perror,'ro');
            %         pause(2)
            %         close all
        end
        if datasum.burstspikes(datasum.steadysweep)>0
            datasum.bursting=datasum.burstspikes(datasum.steadysweep);
        else
            datasum.bursting=0;
        end
        datasum.aphwchange=0;
        datasum.aphwchangevsapnum=0;
        datasum.apamplchange=0;
        datasum.apamplchangevsapnum=0;
        
        if datasum.apnum(datasum.steadysweep)-datasum.burstspikes(datasum.steadysweep)>2
            [p] = polyfit(data.HH.(['sweep',num2str(datasum.steadysweep)]).apmaxtime(datasum.burstspikes(datasum.steadysweep)+1:end),data.HH.(['sweep',num2str(datasum.steadysweep)]).halfwidth(datasum.burstspikes(datasum.steadysweep)+1:end),1);
            datasum.aphwchange=p(1);
            [p] = polyfit(1:length(data.HH.(['sweep',num2str(datasum.steadysweep)]).halfwidth(datasum.burstspikes(datasum.steadysweep)+1:end)),data.HH.(['sweep',num2str(datasum.steadysweep)]).halfwidth(datasum.burstspikes(datasum.steadysweep)+1:end),1);
            datasum.aphwchangevsapnum=p(1);
            [p] = polyfit(data.HH.(['sweep',num2str(datasum.steadysweep)]).apmaxtime(datasum.burstspikes(datasum.steadysweep)+1:end),data.HH.(['sweep',num2str(datasum.steadysweep)]).apamplitude(datasum.burstspikes(datasum.steadysweep)+1:end),1);
            datasum.apamplchange=p(1);
            [p] = polyfit(1:length(data.HH.(['sweep',num2str(datasum.steadysweep)]).apamplitude(datasum.burstspikes(datasum.steadysweep)+1:end)),data.HH.(['sweep',num2str(datasum.steadysweep)]).apamplitude(datasum.burstspikes(datasum.steadysweep)+1:end),1);
            datasum.apamplchangevsapnum=p(1);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%FIRING PATTERN!!!!!!!!!!!!
        %%%%%%%%%%%%%%%%%%%%%%%%Ingerelhet�s�g%%%%%%%%%%%%%
        %     datasum.ingerelhetoseg1=datasum.apnum(datasum.reobasesweep)/20;
        %     datasum.ingerelhetoseg2=0;
        %     datasum.ingerelhetoseg3=0;
        [p] = polyfix([0:20:(datasum.steadysweep-datasum.reobasesweep+1)*20],[0,datasum.apnumfromreobase(1:(datasum.steadysweep-datasum.reobasesweep+1))],0,0,1);
        datasum.excitability=p(1);
        %     if datasum.reobasesweep<length(datasum.apnum)
        %         [p] = polyfit([0:20:40],[0,datasum.apnumfromreobase(1:2)],1);
        %         datasum.ingerelhetoseg2=p(1);
        %         if datasum.reobasesweep+1<length(datasum.apnum)
        %             [p] = polyfit([0:20:60],[0,datasum.apnumfromreobase(1:3)],1);
        %             datasum.ingerelhetoseg3=p(1);
        %         end
        %     end
        
        %%%%%%%%%%%%%%%%%%%%%%%%Ingerelhet�s�g%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%RAMP%%%%%%%%%%%%%%%%%%%%
        datasum.ramp=0;
        datasum.rheobaseramp=0;
        if isfield(data.HH, 'ramp') && datasum.reobasesweep>1 && size(data.HH.ramp,1)>datasum.reobasesweep-2
            datasum.ramp=data.HH.ramp(datasum.reobasesweep-1,1);
        end
        if isfield(data.HH, 'ramp') && size(data.HH.ramp,1)>datasum.reobasesweep-1 && data.HH.(['sweep',num2str(datasum.reobasesweep)]).apmaxtime(1)>.3
            datasum.rheobaseramp=data.HH.ramp(datasum.reobasesweep,1);
        end
        if datasum.reobasesweep>1
            datasum.rampnew=data.HH.rampnew(datasum.reobasesweep-1,1);
        else
            datasum.rampnew=0;
        end
        datasum.rheobaserampnew=data.HH.rampnew(datasum.reobasesweep,1);
        %%%%%%%%%%%%%%%%%%%%%%%%%RAMP%%%%%%%%%%%%%%%%%%%%%
        
        
    end
else
    datasum.aphw=NaN;
    datasum.apampl=NaN;
    datasum.threshold=NaN;
    datasum.apwidth=NaN;
    datasum.apstartenddiff=NaN;
    datasum.dvmax=NaN;
    datasum.dvmaxv=NaN;
    datasum.dvmin=NaN;
    datasum.dvminv=NaN;
    datasum.steadyaphws=NaN;
    datasum.steadyapampls=NaN;
    datasum.steadythresholds=NaN;
    datasum.steadyapwidths=NaN;
    datasum.steadyapstartenddiffs=NaN;
    datasum.steadyahpampls=NaN;
    datasum.steadyahpwidths=NaN;
    datasum.steadyadpampls=NaN;
    datasum.steadyadpwidths=NaN;
    datasum.steadyahpslowampls=NaN;
    datasum.steadyahpslowwidths=NaN;
    datasum.aphwchange=NaN;
    datasum.aphwchangevsapnum=NaN;
    datasum.apamplchange=NaN;
    datasum.apamplchangevsapnum=NaN;
    datasum.ingerelhetoseg=NaN;
    
    datasum.steadysweep=NaN;
    datasum.reobasesweep=NaN;
    
end
end