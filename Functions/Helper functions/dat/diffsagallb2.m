function dsag = diffsagallb2(params)

global t vm

ntr=(length(params)-4)/2;

d1 = params(3:ntr+2);
t0 = params(ntr+3);
tau1 = params(1);
tau2 = params(2);
v0 = params(ntr+4);
v1 = params(ntr+5:2*ntr+4);

for i=1:ntr,
    dsag(1,(i-1)*length(t)+1:i*length(t)) = vm(:,i) - double_exp(t,d1(i),t0,tau1,tau2,v0,v1(i));
end
