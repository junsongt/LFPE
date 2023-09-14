function res = NEQ2 (rm, rn, rd, E, N)
%RM is the mass ratio 
%RN is the number density ratio
%E = E*/kt the reduced threshold energy in the line of ce(2*n)ers cross
%section 
%N is the number of Sonine Polynomial expansion retained, note the program
%will not work if N is greater than 4. 



% WARNING: 
% xmat has 1 based index;
% a and b have 1 based index;
% a11, a12, a21, a22 have 1 based index;
% alpha has 1 based index;
% rhs has 1 based index;
% aa has 0 based index; 

xmat = zeros(2*N, 2*N);

xn1 = rn;
xn2 = 1;
d1 = rd;
d2 = 1;
xm1 = rm / (1 + rm);
xm2 = 1 - xm1;


% AA collision  matrix
a = mat11(1/2);
% AB collision  matrix
b = mat12(1/2);
% Self collision  matrix
c = 2 * (a + b);


%Get the 4 AB collision operator matrix eleme(2*n)s 
a11 = mat11(xm1);
a22 = mat11(xm2);
a12 = mat12(xm2);
a21 = mat12(xm1);


%-----------------------
% nn = N + 5;
% nnm1 = nn - 1; 
% a = [];
% alpha = zeros;
% %alpha2 = []; 
% s(1,1)=1;
% M1=1/(rm+1);
% 
% M2=rm/(1+rm);
% n = N;
%-----------------------
n = N; 


% Set up the main matrix
% Divided eqs. by Q12*N1*N2; Factor in fro(2*n) of A(I,J)
% is N1*N1*Q/(N1*N2*Q12)

% Set up upper two blocks of the matrix
for i = 2:n
    for j = 1:n
        fac1 = (xn1 / xn2) * sqrt(xm2 / 2) * (2 * d1 / (d1 + d2))^2;
        xmat(i, j) = fac1 * c(i, j) + a11(i, j);
    end
    for j = (n+1):(2*n)
        xmat(i, j) = a12(i, j - n);
    end
end

% Set up lower two blocks of the matrix
for i = (n+1):(2*n)
    for j = 1:n
        xmat(i, j) = a21(i - n, j);
    end
    for j = (n+1):(2*n)
        fac2 = (xn2 / xn1) * sqrt(xm1 / 2) * (2 * d2 / (d1 + d2))^2;
        xmat(i, j) = fac2 * c(i - n, j - n) + a22(i - n, j - n);
    end
end

% Now define 1st equation by condition on a1 coefficie(2*n)s
% xmat(1, :) = zeros(1, (2*n));
% xmat(1, 1) = xn1;
% xmat(1, (n+1)) = xn2;

for i = 1:(2*n)
    xmat(1,i) = 0;
end
xmat(1, 1) = xn1;
xmat(1, (n+1)) = xn2;


% Set up the inhomogeneous column vector
% [alpha, aa] = alf(alpha, aa);
result = alf2(E, N, xm1, xm2, xn1, xn2);
alpha = result.alpha;
A = result.A;

for i = 1:n
   rhs(i) = alpha(1,i);
end
for i = (n+1):(2*n)
   rhs(i) = alpha(2,i-n);
end
rhs(1) = 0;



% Solve the set of simultaneous equations
a = linsolve(xmat, rhs');


% Calculate eta
sum = 0;
for i = 1:n
    sum = sum + a(i) * xm2^i * A(i+1) / A(1); % index!!!
end

for i = (n+1):(2*n)
    sum = sum + a(i) * xm1^(i-n) * A(i-n+1) / A(1); % index!!!
end

% disp([' Nonequilibrium correction, eta = ', num2str(-sum)]);
% % t1 = 1 + rhs(1);
% % t2 = 1 + rhs((n+1));
% % disp([' Temperatures of A and B = ', num2str(t1), ' ', num2str(t2)]);

res.eta = -sum;
res.xmat = xmat;
res.rhs = rhs;
res.A = A;
res.a = a';
res.xm1 = xm1;
res.xm2 = xm2;
res.xn1 = xn1;
res.xn2 = xn2;

end






