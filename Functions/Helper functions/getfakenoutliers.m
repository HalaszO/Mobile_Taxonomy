function outID=getfakenoutliers(datain,values)%, iv, datasum, data)
actualdata=cell2mat(datain(2:end,2:end));
ID=datain(2:end,1);
header=datain(1,2:end);
outs=0;
for i=1:length(values)%size(actualdata,2)
    [new,idx,value]=deleteoutliers(actualdata(:,find(strcmp(header,char(values(i))))));
    for j=1:length(idx)
        outs=outs+1;
        outID(outs,1)={[char(ID(idx(j)))]};
        outID(outs,2)=values(i);
        outID{outs,3}=value(j);
%         plotiv(['elfiz',char(ID(idx(j)))], iv, datasum, data, 1, 1000)    
%         h = gcf;
%         set(h,'name',[char(datain(1,i+1)),' actual:',num2str(value(j)),'  average:' ,num2str(mean(actualdata(:,i)))],'numbertitle','off')
%         figure;
%         subplot(1,2,1);
%         hist(actualdata(:,i))
%         subplot(1,2,2);
%         hist(new)
%         h = gcf;
%         set(h,'name',[char(datain(1,i+1)),' actual:',num2str(value(j)),'  average:' ,num2str(mean(actualdata(:,i)))],'numbertitle','off')
%         pause
%         close all;
    end
end
end