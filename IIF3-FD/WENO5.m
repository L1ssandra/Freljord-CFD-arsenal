function h = WENO5(a,b,c,d,e)

epsilon = 1e-5;

h0 = (2*a - 7*b + 11*c)/6;
h1 = (-b + 5*c + 2*d)/6;
h2 = (2*c + 5*d - e)/6;

beta0 = 13/12*(a - 2*b + c)^2 + 1/4*(a - 4*b + 3*c)^2;
beta1 = 13/12*(b - 2*c + d)^2 + 1/4*(b - d)^2;
beta2 = 13/12*(c - 2*d + e)^2 + 1/4*(3*c - 4*d + e)^2;

w0 = 1/(epsilon + beta0)^2;
w1 = 6/(epsilon + beta1)^2;
w2 = 3/(epsilon + beta2)^2;
S = w0 + w1 + w2;
w0 = w0/S;
w1 = w1/S;
w2 = w2/S;

h = w0*h0 + w1*h1 + w2*h2;

end