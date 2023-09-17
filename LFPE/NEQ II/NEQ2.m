% Author: Junsong Tang
% This function is based on Shizgal, 1971(Paper II), which produces an object containing:
% eta
% 4-block main matrix
% alpha quantity
% A integral(moment)
% coefficient vector: a


function res = NEQ2 (rm, rn, rd, E, N)
%rm is the mass ratio 
%rn is the number density ratio
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

% rhs is essentially alpha quantity
for i = 1:n
   rhs(i) = alpha(1,i);
end
for i = (n+1):(2*n)
   rhs(i) = alpha(2,i-n);
end
rhs(1) = 0;



% Solve the set of simultaneous equations to get coeffcients vector: a
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

end






%==========================================================================
% Helper functions

function result = alf2(E, n, xm1, xm2, xn1, xn2)
% WARNING: 
% s has 0 based index;
% xk has 0 based index;
% xj has 0 based index;
% A has 0 based index;
% alpha has 1 based index;


s = zeros(21,21);
xk = zeros(21,1);
xj = zeros(21,1);

nn = n + 5;
nnm1 = nn - 1;


% Calculate the coefficients for the Sonine polynomials s(i,n)
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

% Calculate the J(k) integrals
expe = exp(-E);
xj(1) = expe / 2;
for i = 2:nn+2
    xi = double(i-1); % index!!!
    xj(i) = E^(i-1) * expe / 2 + xi * xj(i-1); % index!!!
    % xj(i) = E^(i) * expe / 2 + i * xj(i-1);
end

% Calculate the K(k) integrals
xk(1) = expe / 2;
for k = 2:nn+1
    xk(k) = xj(k+1) - E * xj(k);
end

% Calculate the A(j) integrals
% a = zeros(nn+1, 1);
for i = 1:nn+1
    sum = 0;
    for k = 1:i+1
        sum = sum + s(i, k) * xk(k);
    end
    A(i) = 4 * sum;
    
end



% Calculate the Alpha(i) quantities
alpha(1, 1) = (xm2 - xn1/(xn1 + xn2)) * A(2); % index!!!
alpha(2, 1) = (xm1 - xn2/(xn1 + xn2)) * A(2); % index!!!
% alpha(1, 1) = (xm2 - xm1*xn1/(xn1 + xn2)) * A(2); % index!!!
% alpha(2, 1) = (xm1 - xm2*xn2/(xn1 + xn2)) * A(2); % index!!!
for k = 2:nn
    alpha(1, k) = xm2^k * A(k+1); % index!!!
    alpha(2, k) = xm1^k * A(k+1); % index!!!
end

result.alpha = alpha;
result.A = A;
result.s = s;
end


% Table II matrix
function a = mat11(xm)
ym = 1 - xm;    
ym2 = ym^2;
ym3 = ym^3;
ym4 = ym^4;
xm2 = xm^2;
xm3 = xm^3;
xm4 = xm^4;
xm5 = xm^5;
xm6 = xm^6;

a(1, 1) = -8 * xm * ym;
a(1, 2) = 4 * xm * ym2;
a(1, 3) = xm * ym3;
% a(1, 4) = xm * ym * (xm3 - 3*xm2 + 3*xm - 1)/ 2; % from paper: a(1, 4) = xm * ym * (xm3 - 3*xm2 + 3*xm - 1)/ 2 = -xm * ym4/ 2
a(1, 4) = xm * ym4/ 2; % from code NEQ2
a(2, 2) = -2 * xm * ym * (15 * xm2 - 18 * xm + 13);
a(2, 3) = xm * ym2 * (35 * xm2 - 30 * xm + 23) / 2;
a(2, 4) = xm * ym * (21 * xm4 - 56 * xm3 + 60 * xm2 - 36 * xm + 11) / 4;
a(3, 3) = -xm * ym * (945 * xm4 - 2100 * xm3 + 2190 * xm2 - 1188 * xm + 433) / 8;
% a(3, 4) = xm * ym * (1155 * xm5 - 3255 * xm4 + 4130 * xm3 - 2970 * xm2 + 1299 * xm - 359) / 16; %from paper
a(3, 4) = -xm * ym * (1155 * xm5 - 3255 * xm4 + 4130 * xm3 - 2970 * xm2 + 1299 * xm - 359) / 16; % from code NEQ2
a(4, 4) = -xm * ym * (15015 * xm6 - 48510 * xm5 + 72135 * xm4 - 62300 * xm3 + 34485 * xm2 - 12102 * xm + 2957) / 32;

for i = 1:3
    for j = i+1:4
        a(j, i) = a(i, j);
    end
end
end



% Table I matrix
function a = mat12(xm)
ym = 1 - xm;
ym2 = ym^2;
ym3 = ym^3;
ym4 = ym^4;
xm2 = xm^2;
xm3 = xm^3;
xm4 = xm^4;

a(1, 1) = 8 * xm * ym;
a(1, 2) = -4 * xm * ym2;
a(1, 3) = -xm * ym3;
a(1, 4) = -xm * ym4 / 2;
a(2, 2) = 30 * xm2 * ym2;
a(2, 3) = -35 * xm2 * ym3 / 2;
a(2, 4) = -21 * xm2 * ym4 / 4;
a(3, 3) = 945 * xm3 * ym3 / 8;
a(3, 4) = -1155 * xm3 * ym4 / 16;
a(4, 4) = 15015 * xm4 * ym4 / 32;

for i = 1:3
    for j = i+1:4
        a(j, i) = a(i, j);
    end
end
end

