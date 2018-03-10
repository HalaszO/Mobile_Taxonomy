function [data] = mHH(iv,dvrs,vrs)
    %BRIDGE BALANCE CHECK - correction did regardless
    %dvrs voltage drop (v0s-vrs)
showprogress = 1;
searchslowevents = 0;
movingaverageforsag = .01;
thresholdvalue = 10;
time_values = iv.time;
sampleinterval = time_values(2)-time_values(1);
%samplingrate=1/(x(2)-x(1)); %sampling rate
xx=abs(time_values-iv.segment(1)/1000);
pulse_start = find(xx==min(xx))+1;
xx = abs(time_values-sum(iv.segment(1:2))/1000);
pulse_end = find(xx==min(xx))-1;

time_values = time_values(pulse_start:pulse_end); % We only care about the values during the current injection

for sweep = 1:iv.sweepnum;

    voltage_values = iv.(['v',num2str(sweep)]);
    voltage_values = voltage_values(pulse_start:pulse_end); % We only care for the values during the current injection
    
   % apmaskfirst = zeros(size(voltage_values));
%     apmaskfirst(find(y+dvrs(sweep)>0))=1;
    apmaskfirst = voltage_values+dvrs(sweep) > 0;  % AP DETECTION THRESHOLD IS 0 mV
    [apmask, apnum] = bwlabelhomemade(apmaskfirst);
    %%%apmask(find(y>0.011))=1;
    %%%[apmask, apnum] = bwlabel(apmask);
    if apnum>800 % inprobable amount of APs 
        apnum=0;
    elseif apnum>0
        %         potentialtrains=[];
        %         for i=1:apnum
        %             if length(y(apmask==i))*sampleinterval>.01
        %                 potentialpeaks=[];
        %                 tempy=y(apmask==i);
        %                 while mean(tempy)>0
        %                     tempeak_h=find(tempy==max(tempy),1,'first');
        %                     tempsearchback_h=tempeak_h;
        %                     while tempsearchback_h>1 && tempy(tempsearchback_h)>tempy(tempsearchback_h-1)
        %                         tempsearchback_h=tempsearchback_h-1;
        %                     end
        %                     tempsearchforward_h=tempeak_h;
        %                     while tempsearchforward_h<length(tempy) && tempy(tempsearchforward_h)>tempy(tempsearchforward_h+1)
        %                         tempsearchforward_h=tempsearchforward_h+1;
        %                     end
        %                     potentialpeakamplitude1=tempy(tempeak_h)-tempy(tempsearchback_h);
        %                     potentialpeakamplitude2=tempy(tempeak_h)-tempy(tempsearchforward_h);
        %                     if potentialpeakamplitude1>.01 && potentialpeakamplitude2>.01
        %                         potentialpeaks(length(potentialpeaks)+1)=tempeak_h;
        %                     end
        %                     tempy(tempsearchback_h:tempsearchforward_h)=0;
        %                 end
        %                 potentialpeaks = sort(potentialpeaks);
        %                 if length(potentialpeaks)>1
        %                     for tempi=2:length(potentialpeaks)
        %                         if (potentialpeaks(tempi)-potentialpeaks(tempi-1))*sampleinterval>.001
        %                             apmaskfirst(find(apmask==i,1,'first')+potentialpeaks(tempi-1)+round((potentialpeaks(tempi)-potentialpeaks(tempi-1))/2))=0;
        %                         end
        %                     end
        %                 end
        %             end
        %         end
        %         [apmask, apnum] = bwlabel(apmaskfirst);
        
        dy = diff(voltage_values)./diff(time_values); % derivative, size: (v_v-1) x 1
        dy0 = dy; 
        dy = mean([0,dy';dy',0])'; % averaging neighbouring elements (first and last values are distorted); size: 1x(dy+1 = v_v)
        yav = mean([voltage_values(2:end)';voltage_values(1:end-1)'])'; % for each dy there's a mean voltage  size: 1x(v_v-1)
        xav = mean([time_values(2:end)';time_values(1:end-1)'])'; % for each dy there's a mean time value
        % NEEDS REVISING FOR ARRAY SIZE
        interspike = 0; %[]
        for ind = 1:apnum
            apmax = max(voltage_values(apmask==ind)); %searching for AP max
            apmax=apmax(1);
            if find((voltage_values == apmax) .* apmask==ind,1,'first')>length(time_values)-50
                apnum=apnum-1; % if the max values are constant, this isn't an AP but an error
            elseif find((voltage_values==apmax) .* apmask==ind,1,'first')<5 && apnum==1
                apnum=0; % Probably detecting initial noise/impedance voltage jumps
            elseif find((voltage_values==apmax).*apmask==ind,1,'first')<5
            % no commands to execute
                
            else
                apmaxtime = time_values((voltage_values==apmax).*apmask==ind);
                apmaxtime = apmaxtime(1);
                apmaxindex = find(time_values==apmaxtime);
                thresh_index = apmaxindex; % searching for threshold starting from the AP apex
                while voltage_values(thresh_index)>0.02 && thresh_index>4 % going under +20 mV 
                    thresh_index=thresh_index-1;
                end
                [~,threshpre] = max(dy(thresh_index-3:apmaxindex-2)); %searches for max derivative in the radius of current th_ind
                thresh_index = threshpre+thresh_index-3;
                while dy(thresh_index)>thresholdvalue && thresh_index>4  %%stepping back 'till we find a threshold in accordance with the value
                    thresh_index=thresh_index-1;
                end
                thresh_index=thresh_index+1;
                threshold=voltage_values(thresh_index); 
                threshtime=time_values(thresh_index);
                % Calculating some correction(?) value (between 0 and 1)
                threshcorrect = 1-abs((thresholdvalue-dy(thresh_index-1))/(dy(thresh_index)-dy(thresh_index-1)));
                data.(['sweep',num2str(sweep)]).threshcorrect(ind) = threshcorrect;
                if data.(['sweep',num2str(sweep)]).threshcorrect(ind) > 1 || data.(['sweep',num2str(sweep)]).threshcorrect(ind)<0
                    data.(['sweep',num2str(sweep)]).threshcorrect(ind) = .5;
                end
                data.(['sweep',num2str(sweep)]).thresh(ind) = threshold-abs(threshcorrect*(voltage_values(thresh_index)-voltage_values(thresh_index-1)));
                data.(['sweep',num2str(sweep)]).threshtime(ind) = threshtime-abs(threshcorrect*(sampleinterval));
                % !!!!!!!!!!!!!!!!!!!!!
%                 if data.(['sweep',num2str(sweep)]).threshtime(sweep)<0
%                     pause;
%                 end 
                apend=find(time_values==apmaxtime)+1; % searching for the end of the AP at the apex
                while voltage_values(apend)>0.01 && apend<length(time_values)-15  %% going under +10mV
                    apend=apend+1;
                end
                while dy(apend)<-thresholdvalue/2 && apend<length(time_values)-15  %% finding -threshold/2
                    apend=apend+1;
                end
                apendd=voltage_values(apend);
                apendtime=time_values(apend);
                apendcorrect = abs((dy(apend)+thresholdvalue/2)/(dy(apend)-dy(apend-1)));
                if isinf(apendcorrect)
                    apendcorrect = .5;
                end
                data.(['sweep',num2str(sweep)]).apendcorrect(ind) = apendcorrect;
                data.(['sweep',num2str(sweep)]).apend(ind)=apendd+abs(apendcorrect*(voltage_values(apend)-voltage_values(apend-1)));
                data.(['sweep',num2str(sweep)]).apendtime(ind)=apendtime-abs(apendcorrect*(sampleinterval));
                
                % Searching for amplitude half-width
                % stepping
                temp=thresh_index;
                while voltage_values(temp)<(apmax-threshold)/2+threshold && temp<length(voltage_values)
                    temp=temp+1;
                end
                aphw_helper = (apmax-threshold)/2+threshold-voltage_values(temp-1);
                data.(['sweep',num2str(sweep)]).aphwstart(ind) = aphw_helper /(voltage_values(temp)-voltage_values(temp-1))*sampleinterval+time_values(temp-1);
                data.(['sweep',num2str(sweep)]).aphwsv(ind)=(apmax-threshold)/2+threshold;
                temp=find(time_values==apmaxtime);
                while voltage_values(temp)>(apmax-threshold)/2+threshold && temp<length(voltage_values)
                    temp=temp+1;
                end
                aphw_helper = (apmax-threshold)/2+threshold-voltage_values(temp-1);
                data.(['sweep',num2str(sweep)]).aphwstop(ind)=aphw_helper/(voltage_values(temp)-voltage_values(temp-1))*sampleinterval+time_values(temp-1);
                data.(['sweep',num2str(sweep)]).apmax(ind)=apmax;
                data.(['sweep',num2str(sweep)]).apmaxtime(ind)=apmaxtime;
                data.(['sweep',num2str(sweep)]).oldthresh(ind)=threshold;
                data.(['sweep',num2str(sweep)]).oldthreshtime(ind)=threshtime;

                data.(['sweep',num2str(sweep)]).oldapend(ind)=apendd;
                data.(['sweep',num2str(sweep)]).oldapendtime(ind)=apendtime;
                data.(['sweep',num2str(sweep)]).apamplitude(ind)=(apmax-threshold)*1000;
                data.(['sweep',num2str(sweep)]).halfwidth(ind)=(data.(['sweep',num2str(sweep)]).aphwstop(ind)-data.(['sweep',num2str(sweep)]).aphwstart(ind))*1000;
                dvmaxlength = find(dy0(thresh_index-3:apend+3)==max(dy0(thresh_index-3:apend+3)),1,'first');
                dvminlength = find(dy0(thresh_index-3:apend+3)==min(dy0(thresh_index-3:apend+3)),1,'first');
                data.(['sweep',num2str(sweep)]).dvmax(ind)=dy0(thresh_index+dvmaxlength-4);
                data.(['sweep',num2str(sweep)]).dvmaxv(ind)=yav(thresh_index+dvmaxlength-4);
                data.(['sweep',num2str(sweep)]).dvmaxt(ind)=xav(thresh_index+dvmaxlength-4);
                data.(['sweep',num2str(sweep)]).dvmin(ind)=dy0(thresh_index+dvminlength-4);
                data.(['sweep',num2str(sweep)]).dvminv(ind)=yav(thresh_index+dvminlength-4);
                data.(['sweep',num2str(sweep)]).dvmint(ind)=xav(thresh_index+dvminlength-4);
                data.(['sweep',num2str(sweep)]).postApexMaxDer(ind)=max(dy0((find(time_values==apmaxtime)+1):apend));

            end
            
        end
    end
    if apnum>0 && ~isfield(data, ['sweep',num2str(sweep)])
        apnum=0;
    end
    if apnum>0
        %%%%%%%%%%%%%%%% DELETING BAD SPIKES
        spikestokill=find(data.(['sweep',num2str(sweep)]).halfwidth<.05 | data.(['sweep',num2str(sweep)]).halfwidth>20 | isnan(data.(['sweep',num2str(sweep)]).halfwidth) | data.(['sweep',num2str(sweep)]).apamplitude<0);
        % HERE
        % If the spike amplitude or halfwidth is not right, the spike is to
        % be discarded
        if length(spikestokill)>100
            spikestokill=1:apnum;
        end
        data.(['sweep',num2str(sweep)]).threshcorrect(spikestokill)=[];
        data.(['sweep',num2str(sweep)]).thresh(spikestokill)=[];
        data.(['sweep',num2str(sweep)]).threshtime(spikestokill)=[];
        data.(['sweep',num2str(sweep)]).apendcorrect(spikestokill)=[];
        data.(['sweep',num2str(sweep)]).apend(spikestokill)=[];
        data.(['sweep',num2str(sweep)]).apendtime(spikestokill)=[];
        data.(['sweep',num2str(sweep)]).aphwstart(spikestokill)=[];
        data.(['sweep',num2str(sweep)]).aphwsv(spikestokill)=[];
        data.(['sweep',num2str(sweep)]).aphwstop(spikestokill)=[];
        data.(['sweep',num2str(sweep)]).apmax(spikestokill)=[];
        data.(['sweep',num2str(sweep)]).apmaxtime(spikestokill)=[];
        data.(['sweep',num2str(sweep)]).oldthresh(spikestokill)=[];
        data.(['sweep',num2str(sweep)]).oldthreshtime(spikestokill)=[];
        data.(['sweep',num2str(sweep)]).oldapend(spikestokill)=[];
        data.(['sweep',num2str(sweep)]).oldapendtime(spikestokill)=[];
        data.(['sweep',num2str(sweep)]).apamplitude(spikestokill)=[];
        data.(['sweep',num2str(sweep)]).halfwidth(spikestokill)=[];
        data.(['sweep',num2str(sweep)]).dvmax(spikestokill)=[];
        data.(['sweep',num2str(sweep)]).dvmaxv(spikestokill)=[];
        data.(['sweep',num2str(sweep)]).dvmaxt(spikestokill)=[];
        data.(['sweep',num2str(sweep)]).dvmin(spikestokill)=[];
        data.(['sweep',num2str(sweep)]).dvminv(spikestokill)=[];
        data.(['sweep',num2str(sweep)]).dvmint(spikestokill)=[];
        data.(['sweep',num2str(sweep)]).postApexMaxDer(spikestokill)=[];
        apnum=apnum-length(spikestokill);
        if apnum>400 % inprobable amount of APs
            apnum=0;
            data.(['sweep',num2str(sweep)]).threshcorrect=[];
            data.(['sweep',num2str(sweep)]).thresh=[];
            data.(['sweep',num2str(sweep)]).threshtime=[];
            data.(['sweep',num2str(sweep)]).apendcorrect=[];
            data.(['sweep',num2str(sweep)]).apend=[];
            data.(['sweep',num2str(sweep)]).apendtime=[];
            data.(['sweep',num2str(sweep)]).aphwstart=[];
            data.(['sweep',num2str(sweep)]).aphwsv=[];
            data.(['sweep',num2str(sweep)]).aphwstop=[];
            data.(['sweep',num2str(sweep)]).apmax=[];
            data.(['sweep',num2str(sweep)]).apmaxtime=[];
            data.(['sweep',num2str(sweep)]).oldthresh=[];
            data.(['sweep',num2str(sweep)]).oldthreshtime=[];
            data.(['sweep',num2str(sweep)]).oldapend=[];
            data.(['sweep',num2str(sweep)]).oldapendtime=[];
            data.(['sweep',num2str(sweep)]).apamplitude=[];
            data.(['sweep',num2str(sweep)]).halfwidth=[];
            data.(['sweep',num2str(sweep)]).dvmax=[];
            data.(['sweep',num2str(sweep)]).dvmaxv=[];
            data.(['sweep',num2str(sweep)]).dvmaxt=[];
            data.(['sweep',num2str(sweep)]).dvmin=[];
            data.(['sweep',num2str(sweep)]).dvminv=[];
            data.(['sweep',num2str(sweep)]).dvmint=[];
            data.(['sweep',num2str(sweep)]).postApexMaxDer=[];
            
        end
        
        if apnum==0
            data=rmfield(data,['sweep',num2str(sweep)]);
        end
        data.NumofDeletedSpikes(sweep) = length(spikestokill);
        
        
        %%%%%%%%%%%%%%%% DELETING BAD SPIKES
        for ind=1:apnum
            if ind>1
                interspike(ind-1)=(data.(['sweep',num2str(sweep)]).apmaxtime(ind)-data.(['sweep',num2str(sweep)]).apmaxtime(ind-1))*1000;
            end
        end
        
        %         if apnum>0
        %%% INTERSPIKE INTERVAL CALCULATION
        if isempty(interspike)==0 
            data.(['sweep',num2str(sweep)]).interspike=[mean(interspike),interspike];
        end
        stepsinoneway=round(.005/sampleinterval); %% Determining the interval for meanspike
        % 5 ms step size or radius
        %%% At 5000Hz (0.0002s sampling), this is 25
        meanspike.(['sweep',num2str(sweep)]).time=([0:sampleinterval:stepsinoneway*sampleinterval])*1000;
        for ind=1:apnum
            meanspike.(['sweep',num2str(sweep)]).(['t',num2str(ind)])=([((-round(stepsinoneway/2)+data.(['sweep',num2str(sweep)]).threshcorrect(ind))*sampleinterval):sampleinterval:stepsinoneway*sampleinterval+data.(['sweep',num2str(sweep)]).threshcorrect(ind)*sampleinterval])*1000;
            % 
            % 25+13+1 num of elements
            % 
            if find(time_values==data.(['sweep',num2str(sweep)]).oldthreshtime(ind))+stepsinoneway>length(time_values)
                meanspike.(['sweep',num2str(sweep)]).(['v',num2str(ind)])=iv.(['v',num2str(sweep)])(-round(stepsinoneway/2)+(find(time_values==data.(['sweep',num2str(sweep)]).oldthreshtime(ind))+pulse_start):pulse_end);
            else
                meanspike.(['sweep',num2str(sweep)]).(['v',num2str(ind)])=iv.(['v',num2str(sweep)])(-round(stepsinoneway/2)+(find(time_values==data.(['sweep',num2str(sweep)]).oldthreshtime(ind)))+pulse_start:(find(time_values==data.(['sweep',num2str(sweep)]).oldthreshtime(ind))+stepsinoneway)+pulse_start);
            end
        end
        %         else
        %             data=rmfield(data,['sweep',num2str(sweep)]);
        %         end
        
        %                 figure;
        %                 for i=1:apnum
        %
        %                     subplot(2,1,1);
        %                     plot(meanspike.(['sweep',num2str(sweep)]).time, meanspike.(['sweep',num2str(sweep)]).(['v',num2str(sweep)]),'-');
        %                     hold on;
        %                     subplot(2,1,2);
        %                     plot(meanspike.(['sweep',num2str(sweep)]).(['t',num2str(sweep)]), meanspike.(['sweep',num2str(sweep)]).(['v',num2str(sweep)]),'-');
        %                     %plot(meanspike.(['sweep',num2str(sweep)]).(['v',num2str(sweep)]) ,[diff(meanspike.(['sweep',num2str(sweep)]).(['v',num2str(sweep)]));0]);
        %                     hold on;
        %
        %
        %                 end
        %                 %xlim([min(meanspike.(['sweep',num2str(sweep)]).time) max(meanspike.(['sweep',num2str(sweep)]).time)]);
        %                 hold off;
        
        
        
        %
        %                 figure;
        %                 plot(x,y);
        %                 hold on
        %                 for i=1:apnum
        %                     plot(data.(['sweep',num2str(sweep)]).apmaxtime(sweep),data.(['sweep',num2str(sweep)]).apmax(sweep),'ro');
        %                     plot(data.(['sweep',num2str(sweep)]).threshtime(sweep),data.(['sweep',num2str(sweep)]).thresh(sweep),'go');
        %                     plot(data.(['sweep',num2str(sweep)]).aphwstart(sweep),data.(['sweep',num2str(sweep)]).aphwsv(sweep),'r.');
        %                     plot(data.(['sweep',num2str(sweep)]).aphwstop(sweep),data.(['sweep',num2str(sweep)]).aphwsv(sweep),'r.');
        %                     plot(data.(['sweep',num2str(sweep)]).apendtime(sweep),data.(['sweep',num2str(sweep)]).apend(sweep),'ko');
        %                     plot( data.(['sweep',num2str(sweep)]).ahptime(sweep),  data.(['sweep',num2str(sweep)]).ahpv(sweep), 'go');
        %                     plot( data.(['sweep',num2str(sweep)]).ahptimeslow(sweep),  data.(['sweep',num2str(sweep)]).ahpvslow(sweep), 'ro');
        %                     plot( data.(['sweep',num2str(sweep)]).adptime(sweep),  data.(['sweep',num2str(sweep)]).adpv(sweep), 'g*');
        %                     plot( data.(['sweep',num2str(sweep)]).adptimeslow(sweep),  data.(['sweep',num2str(sweep)]).adpvslow(sweep), 'r*');
        %                 end
        %                 hold off;
        %figure;
        %plot(apmask);
        %  figure;
        %  plot(data.(['sweep',num2str(sweep)]).apamplitude, data.(['sweep',num2str(sweep)]).halfwidth);
        %  hold on
        %  for q=1:apnum
        %      text(data.(['sweep',num2str(sweep)]).apamplitude(q), data.(['sweep',num2str(sweep)]).halfwidth(q),num2str(q));
        %  end
        % hold off
    end
    
    data.apnum(sweep)=apnum;
    if iv.current(sweep)>0 && apnum==0 % Searching for hump and ramp with 10 ms moving average
        if ~exist('sag','var');
            sag(1,:)=moving(time_values(1:length(find(time_values<0.4+iv.segment(1)/1000))),round(movingaverageforsag/sampleinterval),'mean');
        end
        sag(2,:)=moving(voltage_values(1:length(find(time_values<0.4+iv.segment(1)/1000))),round(movingaverageforsag/sampleinterval),'mean');
        [data.vhump(sweep),temp]=max(sag(2,:));
        data.thump(sweep)=time_values(temp);
        rampstart=find(time_values<.1+iv.segment(1)/1000 ,1,'last');
        data.ramp(sweep,1:2) = polyfit(sag(1,rampstart:end),sag(2,rampstart:end),1);
        %         subplot(1,3,3);
        %         plot(sag(1,:), sag(2,:));
        %         hold on
        %         plot(x, y);
        %         plot(x(find(sag(2,:)==data.vhump(sweep),1,'first')),data.vhump(sweep),'ro');
    elseif  iv.current(sweep)>0 && apnum>0 && sum(sum([data.(['sweep',num2str(sweep)]).apmaxtime<time_values(pulse_start)+.1; data.(['sweep',num2str(sweep)]).apmaxtime>time_values(pulse_start)+.2]))==apnum
        if sum(data.(['sweep',num2str(sweep)]).apmaxtime>.2)>0 && sum(data.(['sweep',num2str(sweep)]).threshtime>.2)>0 && sum((data.(['sweep',num2str(sweep)]).apmaxtime<.2).*(data.(['sweep',num2str(sweep)]).apmaxtime>.1))==0
            rampstart=round(.1/sampleinterval);
            rampend=find(data.(['sweep',num2str(sweep)]).threshtime(find(data.(['sweep',num2str(sweep)]).threshtime>.1,1,'first'))<time_values,1,'first')-round(.01/sampleinterval);
        else
            rampstart=round(.1/sampleinterval);
            rampend=length(time_values);
        end
        data.ramp(sweep,1:2)=polyfit( time_values(rampstart:rampend),voltage_values(rampstart:rampend),1);
    end
    rampstart=round(.2/sampleinterval);
    data.rampnew(sweep,1:2)=polyfit( time_values(rampstart:end),voltage_values(rampstart:end),1);
       
    if showprogress==1
        progressbar([],[],[],sweep/iv.sweepnum,[]);
    end
end
if exist('meanspike','var')
    data.meanspike=meanspike;
end

if ~isfield(data, 'ramp')
    data.ramp=[0,0];
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%AHP �S ADP KERS�S
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%begin, de nem
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%mindenhol
if any(data.apnum>0) && iv.current(end)>=0
    reobasesweep=find(data.apnum>0,1);
%     if isempty(reobasesweep)
%         reobasesweep=2;
%     end
%     steadysweep=find(abs(data.apnum-8)==min(abs(data.apnum-8)),1,'first');
%     if data.apnum(reobasesweep)>8
%         steadysweep=reobasesweep;
    if any(data.apnum>8)
        steadysweep=find(data.apnum>8,1);
    else
        [~,steadysweep]=max(data.apnum);
    end
%     while data.apnum(steadysweep)<8 & length(data.apnum)>steadysweep %& datasum.apnum(datasum.steadysweep)<16
%         steadysweep=steadysweep+1;
%     end
else
    reobasesweep=0;
    steadysweep=0;
end
data.reobasesweep=reobasesweep;
data.steadysweep=steadysweep;

%[~,maxapsweep]=max(data.apnum);
for sweep=1:iv.sweepnum;
    if sweep==reobasesweep || sweep==steadysweep || and(sweep==iv.sweepnum, data.apnum(end)>0)
        voltage_values=iv.(['v',num2str(sweep)]);
        voltage_values=voltage_values(pulse_start:pulse_end);
        apnum=data.apnum(sweep);
        yahp=moving(voltage_values,round(.001/sampleinterval),'mean');     % 1ms moving average
        if searchslowevents==1;
            yahpslow=moving(voltage_values,round(.005/sampleinterval),'mean');     %5ms m. a.
        end
        
        stepsize=round(.002/sampleinterval); % 2 ms step size
        for ind=1:apnum % Analyzing interspike events (AHP&ADP)
            if ind==apnum
                max_length=length(time_values)-stepsize*10;%find(abs(x-sum(iv.segment(1:2))/1000)==min(abs(x-sum(iv.segment(1:2))/1000)))-50;
            else
                max_length=find(abs(time_values-data.(['sweep',num2str(sweep)]).threshtime(ind+1))==min(abs(time_values-data.(['sweep',num2str(sweep)]).threshtime(ind+1))),1,'first')-stepsize*2;
                % Going until we're 2*stepsize close to the next threshold point
                while max_length<5   %maxlength should be at least 5 
                    max_length=max_length+stepsize;
                end
            end
            
            %             if find(x==data.(['sweep',num2str(sweep)]).apmaxtime(sweep))>maxhossz-50 %%rebound spike eset�n
            %                 maxhossz=length(x)-50;
            %             end
            ahp=find(abs(time_values-data.(['sweep',num2str(sweep)]).apendtime(ind))==min(abs(time_values-data.(['sweep',num2str(sweep)]).apendtime(ind))),1,'first');
            % The index of apend
            if ahp>max_length
                ahp=max_length;
            end
            ahpp=ahp;
            while yahp(ahp+stepsize)<yahp(ahp) && ahp<max_length
                ahp=ahp+stepsize;
            end
            ahp=ahp+stepsize;
            ahp=ahpp + find(voltage_values(ahpp:ahp)==min(voltage_values(ahpp:ahp)) ,1,'first');
            adp=ahp; % The search for ADP starts at AHP
            while (yahp(adp+stepsize)>yahp(adp) || yahp(adp+2*stepsize)>yahp(adp) || yahp(adp+3*stepsize)>yahp(adp)) && adp<max_length-stepsize  %x(adp)< data.(['sweep',num2str(sweep)]).threshtime(i+1)
                adp=adp+stepsize;
            end
            % Seems a bit brute force/arbitrary
            adp=adp+stepsize;
            adp=ahp + find(voltage_values(ahp:adp)==max(voltage_values(ahp:adp)) ,1,'first');
            
            %%
            if adp<max_length && searchslowevents==1; %Searching for slow components
                [~,ahpslow]=min(yahp(adp:max_length));
                clear lol;
                ahpslow=ahpslow+adp;
                
                if not(ind==apnum)
                    [~,adp]=max(voltage_values(ahp:ahpslow));
                    clear lol;
                    adp=ahp+adp;
                end
                
                adpslow=ahpslow+2;
                while adpslow+stepsize*10<max_length && yahpslow(adpslow+stepsize*5)>yahpslow(adpslow)
                    adpslow=adpslow+stepsize*5;
                end
                
                adpslow=ahpslow + find(voltage_values(ahpslow:adpslow)==max(voltage_values(ahpslow:adpslow)) ,1,'first');
                
                if ahp>max_length
                    ahp=max_length;
                end
                if adp>max_length
                    adp=max_length;
                end
                if ahpslow>max_length
                    ahpslow=max_length;
                end
                if adpslow>max_length
                    adpslow=max_length;
                end
                data.(['sweep',num2str(sweep)]).ahptimeslow(ind)=time_values(ahpslow);
                data.(['sweep',num2str(sweep)]).ahpvslow(ind)=voltage_values(ahpslow);
                data.(['sweep',num2str(sweep)]).adptimeslow(ind)=time_values(adpslow);
                data.(['sweep',num2str(sweep)]).adpvslow(ind)=voltage_values(adpslow);
            else
                data.(['sweep',num2str(sweep)]).ahptimeslow(ind)=NaN;
                data.(['sweep',num2str(sweep)]).ahpvslow(ind)=NaN;
                data.(['sweep',num2str(sweep)]).adptimeslow(ind)=NaN;
                data.(['sweep',num2str(sweep)]).adpvslow(ind)=NaN;
            end
            %%
            data.(['sweep',num2str(sweep)]).ahptime(ind)=time_values(ahp);
            data.(['sweep',num2str(sweep)]).ahpv(ind)=voltage_values(ahp);
            data.(['sweep',num2str(sweep)]).adptime(ind)=time_values(adp);
            data.(['sweep',num2str(sweep)]).adpv(ind)=voltage_values(adp);
            data.(['sweep',num2str(sweep)]).maxtime(ind)=time_values(max_length);
            
        end
        

    end
    
     %%%% OFFSET VOLTAGE PART
    if data.apnum(sweep)>0
        data.(['sweep',num2str(sweep)]).apmax_corrected=data.(['sweep',num2str(sweep)]).apmax + dvrs(sweep);
        data.(['sweep',num2str(sweep)]).thresh_corrected=data.(['sweep',num2str(sweep)]).thresh + dvrs(sweep);
        
        data.(['sweep',num2str(sweep)]).dvmaxv_corrected=data.(['sweep',num2str(sweep)]).dvmaxv + dvrs(sweep);
        data.(['sweep',num2str(sweep)]).dvminv_corrected=data.(['sweep',num2str(sweep)]).dvminv + dvrs(sweep);
    end
    if iv.current(sweep)>0 && data.apnum(sweep)==0 % calculating hump
        dvhump=vrs(sweep)-data.vhump(sweep);
        data.rhump(sweep)=-dvhump/(iv.current(sweep))*1000000;
    end
    
    if showprogress==1
        progressbar([],[],[],[],sweep/iv.sweepnum);
    end

% Sag extraction
%sagstats = sag_extractor(iv, sweep);    
    
    
end







end

% iterspikesum=0;
% halfwidthsum=0;
% %figure;
% for sweep=1:iv.sweepnum; %rajzolgatom a spikekat
%     if data.apnum(sweep)>0 %els? spikeket �tlagolom
%         % plot(meanspike.(['sweep',num2str(sweep)]).(['v',num2str(1)]),'-');
%         %  hold on;
%         if length(data.(['sweep',num2str(sweep)]).interspike)>2
%             iterspikesum=[iterspikesum,data.(['sweep',num2str(sweep)]).interspike(2:end)];
%             halfwidthsum=[halfwidthsum,data.(['sweep',num2str(sweep)]).halfwidth(2:end)];
%
%
%         end
%     end
% end
% hold off;
% %   interspikesum(1)=[];
% %    halfiwdthsum(1)=[];
% % figure;
% %plot(iterspikesum(2:end),halfwidthsum(2:end));


