function f = func(x,K)
psiFeature = [sin(x(1)); cos(x(1)); 1];
phiFeature = [sin(x(2)); cos(x(2)); 1];
f = psiFeature'*K*phiFeature;
end
