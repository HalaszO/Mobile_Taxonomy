function [out, datasum]=getdatafromdatasum_sz(fnames, datasum, valuenames)
%out=getdatafromdatasumr(fnames, datasum, data, iv,
%{'SAG','HUMP','Rin','v0','reobase','taumin','p','ingerelhetoseg','ramp','aphw','apampl','steadyaphw','steadyapampl','steadytreshold', 'steadyapwidth','apstartenddiff','steadyahpampl','steadyahpwidth','steadyadpampl','steadyadpwidth','steadyahpslowampl','steadyahpslowwidth' });
% make table from structure

for j=1:length(valuenames)
    ALL=0;
    valuename=char(valuenames(j));
    for i=1:length(fnames)
        nev=fnames{i};
        nev=['elfiz',nev(1:end-4)];

        ALL=ALL+1;
        ALLID{ALL}=nev(6:end);
        if isfield(datasum.(nev),valuename)
            disp(nev)
            disp(valuename)
            disp(datasum.(nev).(valuename))
            ALLdata(ALL,j)=datasum.(nev).(valuename);
        else
            ALLdata(ALL,j)=NaN;
        end
    end
end

out=['ID',valuenames;ALLID',num2cell(ALLdata)];
