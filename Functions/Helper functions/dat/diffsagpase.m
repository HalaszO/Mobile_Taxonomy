function dsag = diffsagpase(params)

global t2 vm2

tau = params(1);
amp = params(2:6);
offset = params(7:11);

for i=1:5,
    dsag(1,(i-1)*length(t2)+1:i*length(t2)) = vm2(:,i) - exp1(t2,amp(i),tau,offset(i));
end
