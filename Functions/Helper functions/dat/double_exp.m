function v = double_exp(t,d1,t0,tau1,tau2,v0,v1)

v = (d1.*exp(-(t-t0)./tau1)+(v0-v1-d1).*exp(-(t-t0)./tau2)+v1).*thr(t0,t)+v0.*thr(-t0,-t);
