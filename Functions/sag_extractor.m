function sagstats = sag_extractor(iv)


% Extraction of cellular parameters from step protocol
% Assumes steps starting at 200 ms, ending at 1000 ms, step current from
% -10 to -100 pA in -10 pA steps, and 5 kHz sampling rate
%
% Output: (1) cellname.txt - list of files processed
%         (2) celldata.dat - 5 cell parameters (input resistance, membrane
%         time constant, estimated capacitance, relative sag amplitude,
%         delay of sag peak) per cell
% 
% Requires the Optimizaton Toolbox

global t vm t2 vm2
time_step = iv.time(2)-iv.time(1); % sampling rate (time step)
sampling_freq = 1/time_step; % sampling frequency (in Hz)
cstep_all = iv.current;  % sequence of current step amplitudes (in pA)
step_start = 1e-3*iv.segment(1);
step_end = 1e-3*sum(iv.segment(1:2));
t=iv.time(1:step_end*sampling_freq);

%    vm=1e-3*datin(1:step_end*sampling_freq,2:end);   
     negstep_ind=(cstep_all<0);
%    ntrall=size(vm,2);
%    negstep_ind=negstep_ind(1:ntrall);
    cstep=cstep_all(negstep_ind);
    ntr=sum(negstep_ind);
    vmatrix = zeros(size(iv.time,1),ntr);
    
for sweep_ind=1:ntr
    vmatrix(:, sweep_ind) = iv.(['v',num2str(sweep_ind)]);
end

    vm=vmatrix(1:step_end*sampling_freq,:);
    
    
    % First approximate double exponential fit with common t0, v0, 2 taus
    
    par0 = [0.06 0.01 1e-3*cstep 0.1 -0.06 -0.06+3e-4*cstep];
    lb = [0.1 0.001 -0.3*ones(1,ntr) 0.05 -Inf*ones(1,ntr+1)];
    ub = [1 0.1 zeros(1,ntr) 0.2 Inf*ones(1,ntr+1)];
    options = optimset('Display','iter','TolFun',1e-8,'MaxIter',30,'MaxFunEvals',1e6);

    [bestpars,~,~,~,~] = lsqnonlin(@diffsagallb2,par0,lb,ub,options);



    d1 = bestpars(3:ntr+2);
    t0 = bestpars(ntr+3);
    tau1 = bestpars(1);
    tau2 = bestpars(2);
    v0 = bestpars(ntr+4);
    v1 = bestpars(ntr+5:2*ntr+4);

 %   bestparsall=[tau1.*ones(ntr,1) tau2.*ones(ntr,1) d1' t0.*ones(ntr,1) v0.*ones(ntr,1) v1'];

   tsagmax = t0 + tau1 .* tau2 ./ (tau1 - tau2) .* log(-(v0-v1-d1)./d1.*tau1./tau2);
     sagmin = double_exp(tsagmax,d1,t0,tau1,tau2,v0,v1);
%    vssamp=v1-v0;
%    vtramp=sagmin-v0;
    
    % More accurate double exponential fit to extract Rin and sag parameters
    
    par0 = [tau1 tau2 d1 t0 v0*ones(1,ntr) min(0,v1-v0*ones(1,ntr))];
    lb = [0.1 0.001 -0.3*ones(1,ntr) 0.05 -Inf*ones(1,2*ntr)];
    ub = [1 0.1 zeros(1,ntr) 0.2 zeros(1,2*ntr)];
    options = optimset('Display','iter','TolFun',1e-8,'MaxIter',1e3,'MaxFunEvals',1e6);

    [bestparsc,~,~,~,~] = lsqnonlin(@diffsagallc3,par0,lb,ub,options);

    d1 = bestparsc(3:ntr+2);
    t0 = bestparsc(ntr+3);
    tau1 = bestparsc(1);
    tau2 = bestparsc(2);
    v0 = bestparsc(ntr+4:2*ntr+3);
    vssampc = bestparsc(2*ntr+4:3*ntr+3);
    v1 = v0 + vssampc;

%     bestparsallc=[tau1.*ones(ntr,1) tau2.*ones(ntr,1) d1' t0.*ones(ntr,1) v0' v1'];

    tsagmaxc = t0 + tau1 .* tau2 ./ (tau1 - tau2) .* log(-(v0-v1-d1)./d1.*tau1./tau2);
    sagminc = double_exp(tsagmaxc,d1,t0,tau1,tau2,v0,v1);
    vtrampc=sagminc-v0;
    relsagampc = vtrampc ./ vssampc - 1;
    sagdelayc = tsagmaxc - 0.1;
    rinc=1e-6*cstep'\vssampc';
    
    relsagamp_med=median(relsagampc(6:end));
    sagdelay_med=1000*median(sagdelayc(6:end));
    
    
    % Approximate single exponential fit to early part of step response
    
    t = t(1:round((step_start+0.032)*sampling_freq));
    vm = vm(1:round((step_start+0.032)*sampling_freq),1:5);
    
    par0 = [0.01 0.1 -0.065 -0.06+3e-4*cstep(1:5)];
    lb = [0 0 -Inf*ones(1,6)];
    ub = [1 Inf*ones(1,7)];
    options = optimset('Display','iter','TolFun',1e-8,'MaxIter',1e3,'MaxFunEvals',1e6);

    [bestpars2,~,residual,~,~] = lsqnonlin(@diffsagpasb2,par0,lb,ub,options);
%
%     for i=1:5,
%         r2=sum(residual((i-1)*length(t)+1:i*length(t)).^2);
%         r2m=sum((vm(:,i)-mean(vm(:,i))).^2);
%         resexp(cell,i)=1-r2./r2m;
%     end
    
 %   t0 = bestpars2(2);
    tau = bestpars2(1);
    v0 = bestpars2(3);
    v1 = bestpars2(4:8);

    % More accurate single expontial fit to determine tau_membrane and C
    
    t2=t(round((step_start+0.002)*sampling_freq):end)-t(round(step_start*sampling_freq)+1);
    vm2=vm(round((step_start+0.002)*sampling_freq):end,:);
    
    lb = [0.003 -Inf*ones(1,5) -Inf*ones(1,5)];
    ub = [0.3 Inf*ones(1,10)];
    par0 = [tau v0*ones(1,5)-v1 v1];
    [bestparse,~,~,~,~] = lsqnonlin(@diffsagpase,par0,lb,ub,options);


%allbestpars=cat(1,bestparsall{:});
%allbestparsc=cat(1,bestparsallc{:});
allbestparse=cat(3,bestparse);

alltaus2=1000*squeeze(allbestparse(1,1,:));
allcs3=alltaus2./rinc';

sagstats.rin = rinc;
sagstats.tau = alltaus2';
sagstats.capacity = allcs3';
sagstats.vsag = sagminc;
sagstats.relsagamp = relsagampc;
sagstats.sagdelay = sagdelayc;

end
