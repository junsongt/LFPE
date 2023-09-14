nmax = 5;
E = 3;
eps = 3;

for np1=1:nmax+3
    n=np1-1;
    for ip1=1:np1
        i=ip1-1;
        son(np1,ip1)=(-1)^i*gamma(n+1.5)/(gamma(i+1.5)*factorial(n-i)*factorial(i));
    end
end

nn = nmax + 5;
nnm1 = nn -1;
s(1, 1) = 1;
for i = 1:nnm1+1
    xi = double(i-1); % index!!!
    s(i+1, 1) = (xi + 1.5) * s(i, 1) / (xi + 1);
end

for i = 2:nn+1
    xi = double(i-1);
    for k = 1:i
        xkk = double(k-1); % index!!!
        s(i, k+1) = -(xi - xkk) * s(i, k) / ((xkk + 1.5) * (xkk + 1));
    end
end




expe = exp(-E);
xj(1) = expe / 2;
for i = 2:nn+2
    xj(i) = E^(i-1) * expe / 2 + (i-1) * xj(i-1); % index!!!
    
end

K(1) = expe / 2;
for k = 2:nn+1
    K(k) = xj(k+1) - E * xj(k);
    
end




es=exp(-eps);
jk(1)=es*(eps+1)/2;
for kp1=2:nmax+3
    jk(kp1)= es*(eps^kp1)/2+kp1*jk(kp1-1);
end

xk(1)=es/2;
for kp1=2:nmax+3
    k=kp1-1;
    xk(kp1)=jk(kp1)-jk(kp1-1)*eps;
end


for i = 1:nn+1
    sum = 0;
    for k = 1:i
        sum = sum + s(i, k) * K(k);
    end
    A(i) = 4 * sum;
    
end


for np1=3:nmax+2
    sa=0;
    for kp1=1:np1
        sa=sa+son(np1,kp1)*xk(kp1);
    end
    a(np1-2)=8*sa/2^(np1-1);
end


% A0 = A(1) = 2*exp(-E) = 4*exp(-E)/2