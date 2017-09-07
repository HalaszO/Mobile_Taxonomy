function  [apmask, apnum] = bwlabelhomemade(apmaskfirst)
% apmaskfirst is a logical array
% apnum returns the number of counted APs while apmask contains the
% evolution of apnum elements as an array
% EXAMPLE
%        [0 1 0 0 1 1 0 1 0 0 0 1 1]
% apnum   0 1 1 1 2 2 2 3 3 3 3 4 4  last value as singleton
% apmask [0 1 1 1 2 2 2 3 3 3 3 4 4] array

apmask=zeros(length(apmaskfirst),1);
apnum=0;
for i=2:length(apmaskfirst)
    if apmaskfirst(i-1)==0 && apmaskfirst(i)==1
        apnum=apnum+1;
        apmask(i)=apnum;
    elseif apmaskfirst(i-1)==1 && apmaskfirst(i)==1
        apmask(i)=apnum;
    end
end
end