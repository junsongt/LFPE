% Author: Junsong Tang
% Description: This function gives an object which contains:
% the eigenvalues and eigenfunctions of Lorentz Fokker-Planck operator with
% reactive collision frequency: R(x);
% 


function result = LFPE_Rx(N, npoly, E, rm)
nmax = N;
% E is epsilon* in reactive cross section
% rm is the mass ratio m1/m2 and
% mass fractions M1 = rm/(1+rm); M2 = 1/(1+rm);

format long e
% -- nmax is the number of quad pts; npoly is to plot P_npoly(x)
% -- Load the alpha_n and beta_n recurrence coefficints
% -- Weight function w(x)=x^2exp(-x^2)

% % If users dont have "ab2-70.dat" file, then uncomment below 
% to call ab_maxwell_p2 to get quadrature points 
% ab = ab_maxwell_p2(nmax,100,100,0,10);
% a = ab(:,1);
% b = ab(:,2);
load ab2-70.dat;
a = ab2_70(1:nmax,1); 
b = ab2_70(1:nmax+2,2);

%==========================================================================
% Part 1: Pseudo-spectral representation of Lorentz Fokker-Planck operator
% (no R(x))
%==========================================================================

rtb = sqrt(b); 
rtb(1) = [];
t = diag(rtb(1:nmax-1),-1) + diag(a) + diag(rtb(1:nmax-1),1);
[f,lambda] = eig(t); 
pt = diag(lambda);
% -- Quadrature points and weights from diagonlization of Jacobi matrix
mu_0 = sqrt(pi)/4; 
wt = mu_0*f(1,:).^2;
v = ones(1,nmax);
w = zeros(1,nmax); 
poly = []; 
polyp = []; 
% Construct the SPEED POLYNOMIAL Basis Set
poly = [poly v./sqrt(mu_0)']';
poly = [poly (pt-a(1))/sqrt(b(2)*mu_0)];
polyp = [polyp w'];
polyp = [polyp (v./sqrt(b(2)*mu_0))'];
% Recurrence relation for the Speed Polynomials
for n=2:nmax-1
    xp = (pt-a(n)).*poly(:,n)/rtb(n)-poly(:,n-1)*rtb(n-1)/rtb(n);
    poly = [poly xp];
    yp = (pt-a(n)).*polyp(:,n)/rtb(n)-polyp(:,n-1)*rtb(n-1)/rtb(n)+...
        poly(:,n)/rtb(n);
    polyp = [polyp yp];
end
% Construct the DERIVATIVE MATRIX D:
for i=1:nmax
    for j=1:nmax
        s = sum(poly(j,:).*polyp(i,:)); %sum over P'_n(x_i)P_n(x_j)
        d(i,j) = wt(j)*s;
        %d(i,j)=s;
    end
end


for i=1:nmax
    for j=1:nmax
        s = 0;
        for l=1:nmax
            % s = s + wt(l)*d(l,i)*d(l,j); % this line will give Rayleigh
            s = s + pt(l)*wt(l)*d(l,i)*d(l,j); 
        end
        L_ps(i,j) = s/sqrt(wt(i)*wt(j));
    end
end


% This gives the well-known eigenvalues of the LFPE (without R(x))
[eigfcns,eigen]=eig(L_ps);
eig1_LFP = diag(eigen)';
eig2_LFP = sort(eig1_LFP,'ascend')';
for i=1:nmax
    if eig1_LFP(i) == eig2_LFP(npoly)
        indxeigf=i;
    else
    end
end





%==========================================================================
% Part 2: Pseudo-spectral Reactive collision frequency R(x) from
% Can. J. Phys. paper to the matrix rep of the LFPE - D'.D
%==========================================================================

% Get the I_n integrals
% rm is the mass ratio m1/m2;
% mass fractions M1 and M2 defined as usual

m1=rm/(1+rm);
m2=1/(1+rm); 
rtm1=sqrt(m1); 
rtm2=sqrt(m2);
eps=E;


% % If R is in unit of reduced energy: x1 (Shizgal & Fitzpatrick, 1978)
% % Calculate gp: G(sqrt(x))
% rtx=sqrt(pt);
% t0=(sqrt(eps)+rtm2*rtx)/rtm1;
% i0=sqrt(pi)*erfc(t0)/2;
% i1=0.5*exp(-t0.^2);
% % i2=0.5*exp(-t0.^2).*t0+i0;
% i2=0.5*(exp(-t0.^2).*t0 + i0);
% g1=m1*i2;
% g2=2*rtm1*rtm2*rtx.*i1;
% g3=(m2*pt-eps).*i0;
% Gp=rtm1*(g1-g2+g3);
% 
% %Calculate gm: G(-sqrt(x))
% t0=(sqrt(eps)-rtm2*rtx)/rtm1;
% i0=sqrt(pi)*erfc(t0)/2;
% i1=0.5*exp(-t0.^2);
% % i2=0.5*exp(-t0.^2).*t0+i0;
% i2=0.5*(exp(-t0.^2).*t0 + i0);
% g1=m1*i2;
% g2=-2*rtm1*rtm2*rtx.*i1;
% g3=(m2*pt-eps).*i0;
% Gm=rtm1*(g1-g2+g3);
% % rcol=(gm-gp)./(rtm1*rtm2*pt);
% % rcol=2*pt.*exp(-pt).*rcol/sqrt(pi);
% 
% % Now get R(x)
% for i = 1:nmax
%     x = pt(i);
%     FacR = 2*x*exp(-x)/sqrt(pi);
%     R(i) = FacR*(Gm(i)-Gp(i))/(sqrt(m1*m2*x));
% end


% If R is in unit of reduce speed x
% Note that x is reduced speed while the above x1 is reduces energy !!!
% Calculate gp: G(x)
x = pt;
t0=(sqrt(eps)+rtm2*x)/rtm1;
i0=sqrt(pi)*erfc(t0)/2;
i1=0.5*exp(-t0.^2);
i2=0.5*(exp(-t0.^2).*t0 + i0);
g1=m1*i2;
g2=2*rtm1*rtm2*x.*i1;
g3=(m2*x.^2-eps).*i0;
Gp=rtm1*(g1-g2+g3);

%Calculate gm: G(-x)
t0=(sqrt(eps)-rtm2*x)/rtm1;
i0=sqrt(pi)*erfc(t0)/2;
i1=0.5*exp(-t0.^2);
i2=0.5*(exp(-t0.^2).*t0 + i0);
g1=m1*i2;
g2=-2*rtm1*rtm2*x.*i1;
g3=(m2*pt-eps).*i0;
Gm=rtm1*(g1-g2+g3);

% Now get R(x)
for i = 1:nmax
    x = pt(i);
    FacR = 2*(x^2)*exp(-x^2)/sqrt(pi);
    R(i) = FacR*(Gm(i)-Gp(i))/(rtm1*rtm2*x);
end


%{
% plot R(x)
disp('R(x)')
pause
plot(pt,R,'-ok')
disp('Line 166 - Reactive Collision Frequency')
pause
%}


% Now the pseudo-spectral representation of LFP with R(x): D_ps
% Add this to the diagonal elements of the Fokker-Planck operator with 
% the factor (de/dr)^2
for i=1:nmax
    for j=1:nmax
        if i == j
            D_ps(i,i) = L_ps(i,i) - R(i);
        else
            D_ps(i,j) = L_ps(i,j);
        end
    end
end

% calculate the eigenspace
[eigfcns,eigen]=eig(D_ps);
eig1=diag(eigen)';
eig2=sort(eig1,'ascend')';



for i=1:nmax
    fac2=sqrt(pt(i)^2*exp(-pt(i)^2)/sqrt(pi));
    fac1=sqrt(wt(i));
    eigf2(i)=eigfcns(i,npoly);
    eigf1(i)=eigfcns(i,indxeigf);
%   fprintf('%9.4f %12.5e %12.4e %12.5e %12.4e %12.4e\n' ,pt(i),...
%       eigf2(i),eigf1(i),eigf1(i)/eigf2(i),fac1,fac2);
%   pause
end


% disp('Correction-to-Rx-rate')
% pause
% %fac=eigf2(14)/abs(eigf2(14));





%=========================================================================
% Part 3: Spectral Represnetation  (In Maxwell basis: x^2 * exp(-x^2))
%=========================================================================

% (1)
% Spectral representation of L (LFP operator, Lo & Shizgal, 2006)

for i = 1:N
    for j = 1:N
        if j == i + 1
            L_sp(i,j) = 2*(i-1)*sqrt(b(i+1));
        elseif j == i - 1
            L_sp(i,j) = 2*(j-1)*sqrt(b(j+1));
        elseif j == i && j > 1
            temp = 0;
            for k = 0:(j-2)
                temp = temp + a(k+1);
            end
            L_sp(i,j) = 2 * (j-1) * a(j) + 2 * temp;
        else
            L_sp(i,j) = 0;
        end
    end
end


% (2)
% spectral representation of R(x): coeffs gamma_n
Gam = [];
for i = 1:N
    gam_i = 0;
    for j = 1:N
        % x = pt(j);
        % Gam_i = Gam_i + 4*pi*(x^2)*R(j)*wt(j)*poly(j, i);
        gam_i = gam_i + R(j)*wt(j)*poly(j, i);
    end
    Gam = [Gam gam_i];
end

% (3)
% spectral representation of R(x) as multiplicative operator in matrix form
% (Shizgal's book Page 152)
for i = 1:N
    for j = 1:N
        R_sp(i,j) = 0;
        for k = 1:N
            % x = pt(k);
            % R_sp(i,j) = R_sp(i,j) + 4*pi*(x^2)*R(k)*wt(k)*poly(k,i)*poly(k,j);
            R_sp(i,j) = R_sp(i,j) + R(k)*wt(k)*poly(k,i)*poly(k,j);
        end
    end
end

% Now the Lorentz Fokker-Planck operator with R(x): D_sp
D_sp = L_sp - R_sp;



% Calculate the eigenspace 
[eigfcns,eigen]=eig(D_sp);
eig3=diag(eigen)';
eig4=sort(eig3,'ascend')';


% R(x) = sum_n (gam_n * Mn(x)) 
% (for the purpose of doublecheck if the sum equals R(x))
R_poly = [];
for i = 1:N
    R_poly_i = 0;
    for j = 1:N
        R_poly_i = R_poly_i + Gam(j) * poly(i, j);
    end
    R_poly = [R_poly, R_poly_i];
end





%==========================================================================
% Part 4: Forming result object
%==========================================================================

% eigenvalues of LFP without R(x) 
result.eigenvalues = eig2_LFP';

% eigenvalues of operator in ps & sp rep: should be equal
result.eigenvalues_Rx_ps = eig2';
result.eigenvalues_Rx_sp = eig4';

result.eigenfunction_Rx = eigf1;

result.lambda0_R = eig2(1);
result.lambda0 = eig2_LFP(1);

result.quadpt = pt';
result.weight = wt;

% ps rep of R(x): Rx should be equal to Rpoly
result.Rx = R;
result.Rpoly = R_poly;

% sp rep of R(x)
result.Gam = Gam;

% Maxwell basis (column j --> [Mj(x1),...Mj(xn)]^T)
result.poly = poly;

% ps and sp representation of operator (with/without R(x))
result.L_ps = L_ps;
result.L_sp = L_sp;
result.D_ps = D_ps;
result.D_sp = D_sp;

% estimate of eta
result.etahat = eig2(1) / (exp(-E));


end