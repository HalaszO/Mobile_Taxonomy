function v = single_exp(t,t0,tau,v0,v1)

v = ((v0-v1).*exp(-(t-t0)./tau)+v1).*thr(t0,t)+v0.*thr(-t0,-t);
