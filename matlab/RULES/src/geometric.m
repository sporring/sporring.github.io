function y = geometric(x,P)

index = find(x>=0);
p = 1/(1+P);
q = P/(1+P);
y = p*q.^x(index);
