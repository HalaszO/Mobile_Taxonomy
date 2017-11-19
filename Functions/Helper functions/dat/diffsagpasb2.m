function dsag = diffsagpasb(params)

global t vm

t0 = params(2);
tau = params(1);
v0 = params(3);
v1 = params(4:8);

for i=1:5,
    dsag(1,(i-1)*length(t)+1:i*length(t)) = vm(:,i) - single_exp(t,t0,tau,v0,v1(i));
end
