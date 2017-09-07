function [data] = mpassive_sag(iv,cellname)
showprogress=1;
plotresults=0;
lasthypsectocount=.1; % 100 ms before the end of the stimulus, see vhyp

time_values = iv.time;

time_step = time_values(2)-time_values(1); % sampling rate (time step)
sampling_freq = 1/time_step; 
hundredmicsstep=floor(.0001*sampling_freq); % rounding, in case   Here it is floor(0.5) = 0
% time step is 0.0002 (sec), if time step is <= 0.0001, the above becomes
% larger than 0

if hundredmicsstep==0
    hundredmicsstep=1;
end

% bridge balance
if isfield(iv,'bridged') && iv.bridged==1
    time_minus_stim = abs(time_values - iv.segment(1)/1000);
    pulse_start = find(time_minus_stim == min(time_minus_stim)) - hundredmicsstep; %???? Rakerdezni       
    % Meg miert nem: (taustart) pulse_start = find(x == iv.segment(1)/1000)?
    
else
%%%%%%%% Searching for RS jump START
    time_minus_stim = abs(time_values-iv.segment(1)/1000);
    pulse_start = find(time_minus_stim==min(time_minus_stim))+1;  %Ideally 1000 (5*200 ms)
    start_radius = round(.003*sampling_freq); %  It's 15 here (0.003*5000)   3 msec radius around pulse_start

    % Selecting radius around stimulus start (time values)
    if pulse_start - start_radius > 0
        time_values = time_values((pulse_start - start_radius):(pulse_start + start_radius));
    else
        time_values=time_values(1:pulse_start+start_radius);
    end
    time_values=moving(time_values,hundredmicsstep,'mean'); % Moving average

    for sweep=1:4 % Why examine only the first 4 sweeps?
        Voltage_values = iv.(['v',num2str(sweep)]);
    
        % Selecting radius around stimulus start (voltage values)
        if pulse_start-start_radius > 0
            Voltage_values = Voltage_values((pulse_start-start_radius):(pulse_start+start_radius)); %csak a szegmens elej�n kotor�sszon
        else
            Voltage_values = Voltage_values(1:pulse_start+start_radius);
        end
    
        Voltage_values = moving(Voltage_values,hundredmicsstep,'mean');
        deriv = (diff(Voltage_values)./diff(time_values)); % Time derivative of voltage 
        % To be edited to max(abs(deriv))
        temp(sweep) = find(deriv == min(deriv),1,'first'); % The index of first min of derivative
    end
    tauold = pulse_start;
    [pulse_start, deriv_mode] = mode(temp); % The 'real' start of the stimulus is at the mode of the deriv. minimum 

    % Reminder: All happens within the radius of the assumed pulse start
    if deriv_mode == 1
        pulse_start = start_radius; % A bit weird
    end
    if tauold-start_radius > 0
        pulse_start = pulse_start+tauold-start_radius+1+floor(hundredmicsstep);
    else
        pulse_start = pulse_start+1+floor(hundredmicsstep);
    end
%%%%%%%% Searching for RS jump END
end

time_values = iv.time;

for sweep=1:iv.sweepnum; 
    Voltage_values = iv.(['v',num2str(sweep)]);
    v0 = mean(Voltage_values((pulse_start-5*hundredmicsstep):(pulse_start-3*hundredmicsstep))); %!!
    vrs = mean(Voltage_values(pulse_start)); % Why use mean? pulse_start is a single index    %!!
    if iv.current(sweep) == 0; % 0 current???
        fNorm = 1000 / (sampling_freq/2);               %# normalized cutoff frequency
        [b,a] = butter(6, fNorm, 'low');  %# 6th order filter   Transfer function
        yfilt = filtfilt(b, a, Voltage_values);
        temp = mean(Voltage_values);
        data.noiselevel = mean(abs(Voltage_values-temp))*1000;
        data.filterednoiselevel = mean(abs(yfilt-temp))*1000;
    end
    time_minus_stim = abs(time_values - (iv.segment(1)+iv.segment(2))/1000); % instead of start, the end of the stimulus
    vhypend = find(time_minus_stim == min(time_minus_stim)); % finding the end of the stimulus (index)
    vhypstart = find(time_values < time_values(vhypend)-lasthypsectocount,1,'last');
    vhyp=mean(Voltage_values(vhypstart:vhypend));
    
    dvrs=v0-vrs;
    dvin=vrs-vhyp;
    
    % PLOTTING
    if plotresults==1
        figure
        subplot(1,3,1) 
        hold on;
        plot(time_values, Voltage_values,'.');
        plot(time_values((pulse_start-5*hundredmicsstep):(pulse_start-3*hundredmicsstep)), Voltage_values((pulse_start-5*hundredmicsstep):(pulse_start-3*hundredmicsstep)), 'or');
        plot(time_values(pulse_start), Voltage_values(pulse_start), 'og');
        xlim([iv.segment(1)/1000-.003,iv.segment(1)/1000+.003]);
        title([cellname(6:end),' segments:',num2str(iv.segment)]);
    end
    
    if iv.current(sweep) < 0 % Searching for sag using 500Hz filtering
        %Volt_val_copy = Voltage_values;
        %         fNorm = 1000 / (samplingrate/2);               %# normalized cutoff frequency
        %         [b,a] = butter(10, fNorm, 'low');  %# 10th order filter
        %y = filtfilt(b, a, y);
        
        % Smoothing with moving average
        if round(.005/time_step) < pulse_start
            Voltage_values = moving(Voltage_values,round(.005/time_step),'mean'); % Here this is 25 
        else
            Voltage_values = moving(Voltage_values,pulse_start-1,'mean');
        end
        
        startof = find(time_values > sum(iv.segment(1:2))/1000,1,'first'); % first index after the stimulus ends 
        %finishof = start_radius(time_values)-round(iv.segment(3)/10000*sampling_freq); % Index of a singleton????
        finishof = sum(iv.segment(1:3))/1000*sampling_freq;
        placeof = startof + find(Voltage_values(startof:finishof) == max(Voltage_values(startof:finishof)),1,'first')-1;
        

        if placeof+round(.01*sampling_freq) <= Voltage_values(finishof)
            [data.vrebound(sweep),temp] = max(Voltage_values(placeof-round(.01*sampling_freq) : placeof+round(.01*sampling_freq))); 
            %50
        else
            [data.vrebound(sweep),temp] = max(Voltage_values(placeof-round(.01*sampling_freq) : end));
        end
        data.trebound(sweep) = time_values(temp(1)+placeof-round(.01*sampling_freq));
        data.dvrebound(sweep) = data.vrebound(sweep)-vhyp-dvrs;
        %%%%%%%
        %[data.vsag(sweep),temp] = min(Voltage_values(pulse_start:start_radius(find(time_values < 0.4 + iv.segment(1)/1000))));
    
        [data.vsag_old(sweep),temp] = min(Voltage_values(pulse_start:find(time_values < 0.4+iv.segment(1)/1000,1,'last')));
        % Searching for sag from 0 to pulse_start+400 ms

        dvsag = vrs-data.vsag_old(sweep);
        % added
        data.relsag_old(sweep) = dvsag;
        data.rsag(sweep) = -dvsag/(iv.current(sweep))*1000000;
        sag_min_ind = temp(1)+pulse_start;
        

        if plotresults==1
            subplot(1,3,3);
            % plot(sag(1,:), sag(2,:));
            hold on
            plot(time_values, Voltage_values);
            plot(time_values(sag_min_ind),data.vsag_old(sweep),'ro');
        end
        

        

        clear temp;
        ytau = Voltage_values(pulse_start:sag_min_ind); % Values for measuring time constant
        data.tauend(sweep) = time_values(sag_min_ind);
        orary = ytau - (ytau(1)-(1-1/exp(1))*(ytau(1)-ytau(end)));
        [~,temp] = min(abs(orary));
        data.tauold(sweep) = temp(1)/sampling_freq*1000;
        m = find(abs(ytau-(ytau(1)-(1-1/exp(1))*(ytau(1)-ytau(end))))<.0005);
        if isempty(m)
            m=find(abs(ytau-(ytau(1)-(1-1/exp(1))*(ytau(1)-ytau(end))))<.001);
            if isempty(m)
                m=find(abs(ytau-(ytau(1)-(1-1/exp(1))*(ytau(1)-ytau(end))))<.005);
            end
        end
        data.tauoldfail(sweep)=(max(m)-min(m))/sampling_freq*1000;
        
    end
      
    data.v0s(sweep)=v0;
    data.vrs(sweep)=vrs;
    data.vhyp(sweep)=vhyp;
    data.dvrs(sweep)=dvrs;
    data.dvin(sweep)=dvin;
    data.rs(sweep)=-dvrs/(iv.current(sweep))*1000000;
    data.rin_old(sweep)=-dvin/(iv.current(sweep))*1000000;
    
    if showprogress==1
        progressbar([],[],sweep/iv.sweepnum,0,0);
    end
    clear temp;
end
if ~isfield('data','noiselevel')
    Voltage_values=iv.v1; % using only sweep 1 
    Voltage_values=Voltage_values(find(time_values>iv.segment(1)/1000+.1,1,'first'):find(time_values>sum(iv.segment(1:2)/1000),1,'first')-1);
    fNorm = 1000 / (sampling_freq/2);               %# normalized cutoff frequency
    [b,a] = butter(6, fNorm, 'low');  %# 6th order filter
    yfilt = filtfilt(b, a, Voltage_values);
    temp=mean(Voltage_values);
    data.noiselevel=mean(abs(Voltage_values-temp))*1000; % deviation  !!!!!!!!!!!!!!!!!!
    data.filterednoiselevel=mean(abs(yfilt-temp))*1000;
end

data.samplingrate=sampling_freq;
data.taustart=pulse_start;
data.v0=mean(data.v0s);

% Sag extraction
        
sagstats = sag_extractor(iv);
data.rin_new = sagstats.rin;
data.tau_new = sagstats.tau;
data.capacity_new = sagstats.capacity;
data.vsag_new = sagstats.vsag;
data.relsag_new = sagstats.relsagamp;
data.sagdelay_new = sagstats.sagdelay;

end