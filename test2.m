Es = linspace(0,10,100);
rm = 1/1000000;
rn = 1/1000000;
rd = 1;
N = 4;
etas = [];
eta1 = [];
eta2 = [];
for j = 1 : length(Es)
    E = Es(j);
    result = testII(rm,rn,rd,E,N);
    eta = result.eta;
    etas = [etas eta];
    % y = (1/4)*exp(-E)*((1/2+E)^2 + (1/30)*(3/4+2*E-E^2));
    % eta1 = [eta1, y];
    % z = (1/32)*exp(-E)*(E^4 - 2*E^3 + (1/2)*E^2 + (1/2)*E + 1/16);
    % eta2 = [eta2, z];
end
plot(Es, etas);
xlabel('$\varepsilon^* / KT$','Interpreter','latex','fontsize',20);
ylabel('$\eta$','Interpreter','Latex','fontsize',20);
hold on
% plot(Es, eta1);
% hold on 
% plot(Es, eta2)
% legend('Shizgal', "Nowa", "Present");


function result = testII(rm,rn,rd,E,N)

xm1 = rm/(1+rm);
xm2 = 1-xm1;

xn1 = rn;
xn2 = 1;

d1 = rd;
d2 = 1;

n = N;

a11 = mat11(xm1);
a22 = mat11(xm2);
a12 = mat12(xm2);
a21 = mat12(xm1);

% AA collision  matrix
a = mat11(1/2);
% AB collision  matrix
b = mat12(1/2);
% Self collision  matrix
c = 2 * (a + b);

% Set up the main matrix
% Divided eqs. by Q12*N1*N2; Factor in front of A(I,J)
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
for i = 1:(2*n)
    xmat(1,i) = 0;
end
xmat(1, 1) = xn1;
xmat(1, n+1) = xn2;

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

a = linsolve(xmat, rhs');

sum = 0;
for i = 1:n
    sum = sum - a(i) * xm2^i * A(i+1) / A(1);
end
for i = (n+1):(2*n)
    sum = sum - a(i) * xm1^(i-n) * A(i-n+1) / A(1);
end

result.eta = sum;
end




