function res = NEQ3 (rm, E, rn, delta, N)
%RM is the mass ratio: n/M
%RN is the number density ratio: n/n_M
%E = E*/kt the reduced threshold energy in the line of centers cross
%section
%N is the number of Sonine Polynomial expansion retained, note the program
%will not work if N is greater than 10.
% delta is (dR/dE)^2
q = 0; 
% xmat = zeros(20, 20);
% xmat = zeros(2*N, 2*N);
% a = zeros(10, 10);
% b = zeros(10, 10);
% a12 = zeros(10, 10);
a11 = zeros(10, 10);
% a22 = zeros(10, 10);
% a21 = zeros(10, 10);
% alpha = zeros(2, 10);
% rhs = zeros(20, 1);
% aa = zeros(11, 1);
a = zeros(N,N);

% WARNING: 
% xmat has 1 based index;
% a and b have 1 based index;
% a11, a12, a21, a22 have 1 based index;
% alpha has 1 based index;
% rhs has 1 based index;
% aa has 0 based index; 


aa = zeros(N,1);

np1 = N+1; 
nt = 2*N;


xm1 = rm / (1 + rm);
% xm1 = xm1+1*10.^(-4); 
xm2 = 1 - xm1;



% a = mat11(a,half);
% b = mat12(b,half);
% a = mat11(half);
% b = mat12(half);

%for i = 1 : 4
%    for j = 1 : 4
%        a(i,j) = 2 * (a(i,j) + b(i,j));
%    end
%end
% a = 2 * (a + b);




%Get the 4 AB collision operator matrix elements 
% a11 = mat11(a11, xm1);
% a22 = mat11(a22, xm2);
% a12 = mat12(a12, xm2);
% a21 = mat12(a21, xm1);


a11 = mat1(xm1);

% a11 = mat1(1/2);
%a22 = mat11(xm2);
%a12 = mat12(xm2);
%a21 = mat12(xm1);


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
% Divided eqs. by Q12*N1*N2; Factor in front of A(I,J)
% is N1*N1*Q/(N1*N2*Q12)

fac1 = 1 * sqrt(xm2 / 2) * (2 * (1/2))^2;
%fac1 = (xn1 / xn2) * sqrt(xm2 / 2) * (2 * d1 / (d1 + d2))^2;





% Set up upper two blocks of the matrix
%for i = 2:n
%    for j = 1:n
%        fac1 = (xn1 / xn2) * sqrt(xm2 / 2) * (2 * d1 / (d1 + d2))^2;
%        xmat(i, j) = fac1 * a(i, j) + a11(i, j);
%    end
%    for j = np1:nt
%        xmat(i, j) = a12(i, j - n);
%    end
%end

% Set up lower two blocks of the matrix
%for i = np1:nt
%    for j = 1:n
%        xmat(i, j) = a21(i - n, j);
%    end
%    for j = np1:nt
%        fac2 = (xn2 / xn1) * sqrt(xm1 / 2) * (2 * d2 / (d1 + d2))^2;
%        xmat(i, j) = fac2 * a(i - n, j - n) + a22(i - n, j - n);
%    end
%end

% Now define 1st equation by condition on a1 coefficients
% xmat(1, :) = zeros(1, nt);
% xmat(1, 1) = xn1;
% xmat(1, np1) = xn2;

%for i = 1:nt
%    xmat(1,i) = 0;
%end
%xmat(1, 1) = xn1;
%xmat(1, np1) = xn2;

%if q == 1
%    disp(' ');
%    disp(' MATRIX');
%    for i = 1:m
%        disp(xmat(i, :));
%    end
%end

% Set up the inhomogeneous column vector
% [alpha, aa] = alf(alpha, aa);
result = alf3(E, N, xm2);
alpha = result.alpha;
A = result.A;

for i = 1:n
   rhs(i) = alpha(i);
end

%for  i = np1:nt
%   rhs(i) = alpha(2,i-n);
%end
%rhs(1) = 0;


% rhs = zeros(1, nt);
% %disp(alpha(1, :))
% rhs(1:nt) = alpha(1: nt);
% rhs(np1:nt) = alpha(2, 1:n);
% rhs(2) = 0;



%a11
% Solve the set of simultaneous equations
% for i = 1:N 
%     a(i,i) = a11(i,i); 
% end
% 
% for i = 1:N
%     for j = i+1:N
%         a(j,i) = a11(j,i); 
%     end
% end
% 
% for i = 1:N
%     for j = i+1:N
%         a(i,j) = a(j, i);
%     end
% end

mtx = a11(1:N, 1:N);
rhs = rhs * rn * delta;

 
a = linsolve(mtx, rhs');
% a = rn * a;

sum = 0;
for i = 1:n
    % sum = sum - rhs(i) * xm2^i * aa(i+1) / aa(1); % index!!!
    sum = sum - a(i) * A(i+1) / A(1);
end

%for i = np1:nt
%    sum = sum + rhs(i) * xm1^(i-n) * aa(i-n+1) / aa(1); % index!!!
%    disp([i, sum]);
%end

disp([' Nonequilibrium correction, eta = ', num2str(sum)]);
%t1 = 1 + rhs(1);
%t2 = 1 + rhs(np1);
%disp([' Temperatures of A and B = ', num2str(t1), ' ', num2str(t2)]);



% Now use 1 - lambda_0/keq to estimate eta
Elastic(2:N+1, 2:N+1) = mtx;
for i = 1:N+1
    Elastic(1,i) = 0;
    Elastic(i,1) = 0;
end


for i = 1:N+1
    for j = 1:N+1
        if i == 1
            R(i, j) = A(j);
        elseif j == 1
            R(i, j) = A(i);
        else
            R(i,j) = 0; % Rij = sum R(xk)*Si(xk)*Sj(xk)*wt(xk)
        end
    end
end




for i = 1:N+1
    for j = 1:N+1
        Ni = 2 * (gamma(i-1+1.5)/factorial(i-1)) / sqrt(pi);
        Nj = 2 * (gamma(j-1+1.5)/factorial(j-1)) / sqrt(pi);
        fac = 1 / (sqrt(Ni) * sqrt(Nj));
        M(i,j) = fac * (Elastic(i,j) - R(i,j));
    end
end




[eigf,eigv]=eig(M);
eigen1 = diag(eigv);
eigen2 = sort(abs(eigen1),'ascend')';



res.eta = sum;
res.a11 = a11;

res.rhs = rhs;
res.a = a;
res.A = A;
res.alpha = alpha;
res.eigen_values = eigen2';

res.mtx = mtx;
res.Elastic = Elastic;
res.Reaction = R;
res.M = M;
res.etahat = 1 - abs(eigen2(1))/A(1);


end






