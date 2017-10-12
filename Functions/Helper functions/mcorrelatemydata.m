function mcorrelatemydata(datain)
colours={'r','g','b','k','m','c','y'};
colors={'red','green','blue','black','magenta','cyan','yellow'};
fields=fieldnames(datain);
temp=datain.(char(fields(1)));
datainsum=temp(1,:);
for i=1:length(fields)
    temp=datain.(char(fields(i)));
    datanums(i)={cell2mat(temp(2:end,2:end))};
    datainsum=[datainsum;temp(2:end,:)];
end
datanumsum=cell2mat(datainsum(2:end,2:end));
valuenames=temp(1,2:end);

[r,p] = corrcoef(datanumsum);
for k=1:length(fields)
    datanum=cell2mat(datanums(k));
    [rrr{k},ppp{k}]=corrcoef(datanum);
end


for i=1:size(r,1)
    for j=1:size(r,1)
        for k=1:length(fields)
            pact=cell2mat(ppp(k));
            ps(k)=pact(i,j);
        end
        if strcmp(char(valuenames(i)),'ISaac') %&& strcmp(char(valuenames(j)),'sampleinterval')%&& sum([p(i,j),ps]<.05)>0  
            text='';
            figure(323);
            for k=1:length(fields)
                datanum=cell2mat(datanums(k));
                plot(datanum(:,i),datanum(:,j),['.',char(colours(k))]);
                hold on;
                [pp] = polyfit(datanum(:,i),datanum(:,j),1);
                plot(datanum(:,i),pp(1)*datanum(:,i)+pp(2),char(colours(k)));
                ract=cell2mat(rrr(k));
                pact=cell2mat(ppp(k));
                if pact(i,j)<.05
                    format short
                    ract=num2str(ract(i,j));
                    pact=num2str(pact(i,j));
                    text=([text,'{\color{',char(colors(k)),'} ',char(fields(k)),' r:',ract,'  p:',pact,' }']);
                end
            end
            text=([text,'{\color{black} SUM: r:',num2str(r(i,j)),'  p:',num2str(p(i,j)),' }']);
            [pp] = polyfit(datanumsum(:,i),datanumsum(:,j),1);
            plot(datanumsum(:,i),pp(1)*datanumsum(:,i)+pp(2),'k');
            h = gcf;
            set(h,'name',[char(valuenames(i)),' vs ',char(valuenames(j)),'  r:',num2str(r(i,j)),'  p:',num2str(p(i,j))],'numbertitle','off');
            legend(text);
            xlabel(char(valuenames(i)));
            ylabel(char(valuenames(j)));
            pause
            clf;
            %close all;
        end
    end
end

end