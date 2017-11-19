% FITSAGALLO5
% Extraction of cellular parameters from step protocol
% Assumes steps starting at 100 ms, ending at 900 ms, step current from
% -10 to -100 pA in -10 pA steps, and 5 kHz sampling rate
%
% Output: (1) cellname.txt - list of files processed
%         (2) celldata.dat - 5 cell parameters (input resistance, membrane
%         time constant, estimated capacitance, relative sag amplitude,
%         delay of sag peak) per cell
% 
% Requires the Optimizaton Toolbox

clear all
close all

global t vm t2 vm2

file_filter='*.dat';    % files to be included in analysis
sampling_rate=10000;      % sampling frequency (in Hz)
cstep_all=-10:-10:-100;  % sequence of current step amplitudes (in pA)
step_start=0.1;
step_end=0.9;

wd=input('Path to data files (RETURN - current directory): ','s');
if strcmp(wd,'')
    wd='.';
elseif ~exist(wd,'dir')
    disp('No such directory.')
    return
end

files=dir(strcat(wd,'/',file_filter));
ncells = length(files);

fid = fopen(strcat(wd,'/cellnames.txt'),'wt');

for cell=1:ncells,
    fname{cell}=files(cell).name;
    cellname{cell}=fname{cell}(1:length(fname{cell})-4);
    
    fprintf(fid,'%s \n',fname{cell});
    
    datin=importdata(strcat(wd,'/',fname{cell}));
    t=1e-3*datin(1:step_end*sampling_rate,1);
    vm=1e-3*datin(1:step_end*sampling_rate,2:end);
    
    negstep_ind=(cstep_all<0);
    ntrall=size(vm,2);
    negstep_ind=negstep_ind(1:ntrall);
    cstep=cstep_all(negstep_ind);
    vm=vm(:,negstep_ind);
    ntr=sum(negstep_ind);
    
    % First approximate double exponential fit with common t0, v0, 2 taus
    
    par0 = [0.06 0.01 1e-3*cstep 0.1 -0.06 -0.06+3e-4*cstep];
    lb = [0.1 0.001 -0.3*ones(1,ntr) 0.05 -Inf*ones(1,ntr+1)];
    ub = [1 0.1 zeros(1,ntr) 0.15 Inf*ones(1,ntr+1)];
    options = optimset('Display','iter','TolFun',1e-8,'MaxIter',100,'MaxFunEvals',1e6);

    [bestpars{cell},resnorm{cell},residual,exitflag,output] = lsqnonlin(@diffsagallb2,par0,lb,ub,options);

    bestpars{cell}
    resnorm{cell}
    output

    d1 = bestpars{cell}(3:ntr+2);
    t0 = bestpars{cell}(ntr+3);
    tau1 = bestpars{cell}(1);
    tau2 = bestpars{cell}(2);
    v0 = bestpars{cell}(ntr+4);
    v1 = bestpars{cell}(ntr+5:2*ntr+4);

    bestparsall{cell}=[tau1.*ones(ntr,1) tau2.*ones(ntr,1) d1' t0.*ones(ntr,1) v0.*ones(ntr,1) v1'];

    tsagmax = t0 + tau1 .* tau2 ./ (tau1 - tau2) .* log(-(v0-v1-d1)./d1.*tau1./tau2);
    sagmin = double_exp(tsagmax,d1,t0,tau1,tau2,v0,v1);
    vssamp=v1-v0;
    vtramp=sagmin-v0;
    
    % More accurate double exponential fit to extract Rin and sag parameters
    
    par0 = [tau1 tau2 d1 t0 v0*ones(1,ntr) min(0,v1-v0*ones(1,ntr))];
    lb = [0.1 0.001 -0.3*ones(1,ntr) 0.05 -Inf*ones(1,2*ntr)];
    ub = [1 0.1 zeros(1,ntr) 0.15 zeros(1,2*ntr)];
    options = optimset('Display','iter','TolFun',1e-8,'MaxIter',1e3,'MaxFunEvals',1e6);

    [bestparsc{cell},resnormc{cell},residual,exitflag,output] = lsqnonlin(@diffsagallc3,par0,lb,ub,options);

    bestparsc{cell}
    resnormc{cell}
    output

    d1 = bestparsc{cell}(3:ntr+2);
    t0 = bestparsc{cell}(ntr+3);
    tau1 = bestparsc{cell}(1);
    tau2 = bestparsc{cell}(2);
    v0 = bestparsc{cell}(ntr+4:2*ntr+3);
    vssampc = bestparsc{cell}(2*ntr+4:3*ntr+3);
    v1 = v0 + vssampc;

    bestparsallc{cell}=[tau1.*ones(ntr,1) tau2.*ones(ntr,1) d1' t0.*ones(ntr,1) v0' v1'];

    tsagmaxc = t0 + tau1 .* tau2 ./ (tau1 - tau2) .* log(-(v0-v1-d1)./d1.*tau1./tau2);
    sagminc = double_exp(tsagmaxc,d1,t0,tau1,tau2,v0,v1);
    vtrampc=sagminc-v0;
    relsagampc = vtrampc ./ vssampc - 1;
    sagdelayc = tsagmaxc - 0.1;
    rinc(cell)=1e-6*cstep(1:5)'\vssampc(1:5)';
    
    relsagamp_med(cell)=median(relsagampc(6:end));
    sagdelay_med(cell)=1000*median(sagdelayc(6:end));
    
    close all
    fh=figure(1);
    subplot(2,2,2)
    plot(cstep,vssampc,'x')
    hold on
    plot(cstep,vtrampc,'o')
    line([cstep(5) 0],[1e-6*rinc(cell)*cstep(5) 0])
    
    subplot(2,2,1)
    set(gcf,'DefaultAxesColorOrder',[1 0 0;0 1 0])
    plot(t,vm,'.')
    hold on
    for i=1:ntr,
        plot(t,double_exp(t,bestpars{cell}(i+2),bestpars{cell}(ntr+3),bestpars{cell}(1),bestpars{cell}(2),bestpars{cell}(ntr+4),bestpars{cell}(ntr+4+i)),'k','LineWidth',2)
        plot(t,double_exp(t,bestparsc{cell}(i+2),bestparsc{cell}(ntr+3),bestparsc{cell}(1),bestparsc{cell}(2),bestparsc{cell}(ntr+3+i),bestparsc{cell}(ntr+3+i)+bestparsc{cell}(2*ntr+3+i)),'m','LineWidth',2)
    end
    title(cellname{cell})
    
    % Approximate single exponential fit to early part of step response
    
    t = t(1:round((step_start+0.032)*sampling_rate));
    vm = vm(1:round((step_start+0.032)*sampling_rate),1:5);
    
    par0 = [0.01 0.1 -0.065 -0.06+3e-4*cstep(1:5)];
    lb = [0 0 -Inf*ones(1,6)];
    ub = [1 Inf*ones(1,7)];
    options = optimset('Display','iter','TolFun',1e-8,'MaxIter',1e3,'MaxFunEvals',1e6);

    [bestpars2{cell},resnorm2{cell},residual,exitflag,output] = lsqnonlin(@diffsagpasb2,par0,lb,ub,options);

    bestpars2{cell}
    resnorm2{cell}
    output

    for i=1:5,
        r2=sum(residual((i-1)*length(t)+1:i*length(t)).^2);
        r2m=sum((vm(:,i)-mean(vm(:,i))).^2);
        resexp(cell,i)=1-r2./r2m;
    end
    
    t0 = bestpars2{cell}(2);
    tau = bestpars2{cell}(1);
    v0 = bestpars2{cell}(3);
    v1 = bestpars2{cell}(4:8);

    % More accurate single expontial fit to determine tau_membrane and C
    
    t2=t(round((step_start+0.002)*sampling_rate):end)-t(round(step_start*sampling_rate)+1);
    vm2=vm(round((step_start+0.002)*sampling_rate):end,:);
    
    lb = [0.003 -Inf*ones(1,5) -Inf*ones(1,5)];
    ub = [0.3 Inf*ones(1,10)];
    par0 = [tau v0*ones(1,5)-v1 v1];
    [bestparse{cell},resnorme{cell},residual,exitflag,output] = lsqnonlin(@diffsagpase,par0,lb,ub,options);
    
    bestparse{cell}
    resnorme{cell}
    output
    
    subplot(2,2,3)
    set(gcf,'DefaultAxesColorOrder',[1 0 0;0 1 0])
    plot(t(round((step_start-0.01)*sampling_rate):end),vm(round((step_start-0.01)*sampling_rate):end,:),'.')
    hold on
    for i=1:5,
        vi(cell,i)=exp1(t(round(step_start*sampling_rate)+1),bestparse{cell}(1+i),bestparse{cell}(1),bestparse{cell}(6+i));
        plot(t(round(step_start*sampling_rate)+1:end),exp1(t(round(step_start*sampling_rate)+1:end)-t(round(step_start*sampling_rate)+1),bestparse{cell}(1+i),bestparse{cell}(1),bestparse{cell}(6+i)),'k','LineWidth',2)
    end
    
    print(fh,'-djpeg',strcat(wd,'/',cellname{cell}))
end

fclose(fid);

allbestpars=cat(1,bestparsall{:});
allbestparsc=cat(1,bestparsallc{:});
allbestparse=cat(3,bestparse{:});

alltaus2=1000*squeeze(allbestparse(1,1,:));
allcs3=alltaus2./rinc';

celldata_out = [rinc' alltaus2 allcs3 relsagamp_med' sagdelay_med'];
dlmwrite(strcat(wd,'/celldata.dat'), celldata_out, 'delimiter','\t', 'precision',4)
