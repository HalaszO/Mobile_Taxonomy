function dsag = diffsagallc3(params)

global t vm

ntr=(length(params)-3)/3;

d1 = params(3:ntr+2);
t0 = params(ntr+3);
tau1 = params(1);
tau2 = params(2);
v0 = params(ntr+4:2*ntr+3);
vssamp = params(2*ntr+4:3*ntr+3);
v1 = v0 + vssamp;

for i=1:ntr,
    dsag(1,(i-1)*length(t)+1:i*length(t)) = vm(:,i) - double_exp(t,d1(i),t0,tau1,tau2,v0(i),v1(i));
end
