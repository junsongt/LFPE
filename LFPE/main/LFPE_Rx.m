function result = LFPE_Rx(N, npoly, E, rm)
nmax = N;
estar = E;
gam = rm;

% estar is epsilon* in reactive cross section
% gam is the mass ratio m1/m2 - mass fractions M1 and M2
%
%xm1=1/(1+gam); xm2=gam/(1+gam);
% This version adds the reactive collision frequency to the
% program Lorentz-FPE.m: This code works and gives the correct
% eigenvalues
format long e
% -- nmax is the number of quad pts; npoly is to plot P_npoly(x)
output=fopen('eig_fpe.dat','wt');
% -- Load the alpha_n and beta_n recurrence coefficints
% -- Weight function w(x)=x^2exp(-x^2)
% ab = ab_maxwell_p2(nmax,100,100,0,10);
% a = ab(:,1);
% b = ab(:,2);
load ab2-70.dat;
a=ab2_70(1:nmax,1); b=ab2_70(1:nmax+2,2);
rtb=sqrt(b); rtb(1)=[];
t=diag(rtb(1:nmax-1),-1)+diag(a)+diag(rtb(1:nmax-1),1);
[f,lambda]=eig(t); pt=diag(lambda);
% -- Quadrature points and weights from diagonlization of Jacobi matrix
mu_0=sqrt(pi)/4; mu_1=0.5;
wt=mu_0*f(1,:).^2;
% -- First three moments of the weight function
rtpi=sqrt(pi); xn=[0:1:nmax]; 
v=ones(1,nmax);w = zeros(1,nmax); poly=[]; polyp=[]; 
% Construct the SPEED POLYNOMIAL Basis Set
poly=[poly v./sqrt(mu_0)']';
poly=[poly (pt-a(1))/sqrt(b(2)*mu_0)];
polyp=[polyp w'];
polyp=[polyp (v./sqrt(b(2)*mu_0))'];
% Recurrence relation for the Speed Polynomials
for n=2:nmax-1
    xp=(pt-a(n)).*poly(:,n)/rtb(n)-poly(:,n-1)*rtb(n-1)/rtb(n);
    poly=[poly xp];
    yp=(pt-a(n)).*polyp(:,n)/rtb(n)-polyp(:,n-1)*rtb(n-1)/rtb(n)+...
        poly(:,n)/rtb(n);
    polyp=[polyp yp];
end
% Construct the DERIVATIVE MATRIX D:
for i=1:nmax
    for j=1:nmax
        s=sum(poly(j,:).*polyp(i,:)); %sum over P'_n(x_i)P_n(x_j)
        d(i,j)=wt(j)*s;
        %d(i,j)=s;
    end
end


for i=1:nmax
    for j=1:nmax
        s=0;
        for l=1:nmax
            % s = s + wt(l)*d(l,i)*d(l,j); % this will give Rayleigh
            s = s + pt(l)*wt(l)*d(l,i)*d(l,j); 
        end
        L_ps(i,j)=s/sqrt(wt(i)*wt(j));
    end
end


% This gives the well-known eigenvalues of the LFPE
[eigfcns,eigen]=eig(L_ps);
eig1_LFP = diag(eigen)';
eig2_LFP = sort(eig1_LFP,'ascend')';
for i=1:nmax
    if eig1_LFP(i) == eig2_LFP(npoly)
        indxeigf=i;
    else
    end
end
imax=nmax;





% =========================================
% Now add the reactive collision frequency R(x) from
% Can. J. Phys. paper to the matrix rep of the LFPE - D'.D
% =========================================
%
%Get the I_n integrals
% gam is the mass ratio m1/m2 - mass fractions M1 and M2
%
m2=1/(1+gam); m1=gam/(1+gam);
rtm1=sqrt(m1); rtm2=sqrt(m2);
eps=estar;


xm1=m1;
xm2=m2;

% % In reduced energy: x1
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
% % ==========
% % Now get R(x)
% % ==========
% for i = 1:nmax
%     x = pt(i);
%     FacR = 2*x*exp(-x)/sqrt(pi);
%     R(i) = FacR*(Gm(i)-Gp(i))/(sqrt(m1*m2*x));
% end



% x is reduced velocity while x1 is reduces energy
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


% ==========
% Now get R(x)
% ==========
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

%
% Add this to the diagonal elements of the Fokker-Planck operator with 
% the factor (de/dr)^2
%

for i=1:nmax
    for j=1:nmax
        if i == j
            D_ps(i,i) = L_ps(i,i) - R(i);
        else
            D_ps(i,j) = L_ps(i,j);
        end
    end
end


[eigfcns,eigen]=eig(D_ps);
eig1=diag(eigen)';
eig2=sort(eig1,'ascend')';

% [eigfcns,eigen]=eig(D);
% eig1=diag(eigen)';
% eig2=sort(eig1,'ascend')';

%{
% print out modified eigenvalues
i=1;
fprintf('%2i   %12.6e \n' ,i, eig2(i));
for i=2:imax
    fprintf('%2i   %12.6f \n' ,i, eig2(i));
end
disp('Modified Eigenvalues')
%}


%{
% check if the index of wanted eigenfunction stay unchanged even if with
% R(x)
pause
disp(indxeigf)
pause
for i=1:nmax
    if eig1(i)==eig2(npoly)
        mod_indxeigf=i;
    else
    end
end
disp(mod_indxeigf)
pause
%}


for i=1:nmax
    fac2=sqrt(pt(i)^2*exp(-pt(i)^2)/sqrt(pi));
    fac1=sqrt(wt(i));
    eigf2(i)=eigfcns(i,npoly);
    eigf1(i)=eigfcns(i,indxeigf);
%   fprintf('%9.4f %12.5e %12.4e %12.5e %12.4e %12.4e\n' ,pt(i),...
%       eigf2(i),eigf1(i),eigf1(i)/eigf2(i),fac1,fac2);
%   pause
end




% % THIS ETA IS WRONG!!!
% eta=(exp(-eps)+eig2(1))/exp(-eps);


% disp('Correction-to-Rx-rate')
% pause
% %fac=eigf2(14)/abs(eigf2(14));

%{
% plot the npoly-th eigenfunctions
plot(pt,eigf1,'-ok','linewidth',1.2,'markersize',7,'markerfacecolor','k')
axis([0 5 -.5 .5])
set(gca,'FontSize',24)
set(gca,'Ytick',[-.4:.2:.4],'linewidth',1.6)
set(gca,'Xtick',[0:1:5],'linewidth',1.6)
str2=num2str(npoly);
str1={'$n = $'};
str=strcat(str1,str2);
text(3,.3,0,str,'Interpreter','latex','fontsize',26)
xlabel('$x$','Interpreter','latex','fontsize',30)
ylabel('$\psi_n(x)$','Interpreter','Latex','fontsize',30)
%}



%=========================================================================
% Spectral representation of L
% D(1,1) = 0;
% for i = 2:N
%     temp = 0;
%     for k = 0:(i-2)
%         temp = temp + a(k+1);
%     end
%     D(i,i) = (i-1) * a(i-1) + temp;
% end
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



% spectral representation of R(x): coeffs gam_n
Gam = [];
for i = 1:N
    Gam_i = 0;
    for j = 1:N
        x = pt(j);
        % Gam_i = Gam_i + 4*pi*(x^2)*R(j)*wt(j)*poly(j, i);
        Gam_i = Gam_i + R(j)*wt(j)*poly(j, i);
    end
    Gam = [Gam Gam_i];
end


% spectral matrix rep of R(x) in Maxwell space
for i = 1:N
    for j = 1:N
        R_sp(i,j) = 0;
        for k = 1:N
            x = pt(k);
            % R_sp(i,j) = R_sp(i,j) + 4*pi*(x^2)*R(k)*wt(k)*poly(k,i)*poly(k,j);
            R_sp(i,j) = R_sp(i,j) + R(k)*wt(k)*poly(k,i)*poly(k,j);
        end
    end
end

D_sp = L_sp - R_sp;


% for i=1:nmax
%     for j=1:nmax
%         if i == j
%             D_sp(i,i) = L_sp(i,i) - Gam(i);
%         else
%             D_sp(i,j) = L_sp(i,j);
%         end
%     end
% end

[eigfcns,eigen]=eig(D_sp);
eig3=diag(eigen)';
eig4=sort(eig3,'ascend')';


% R(x) value
R_poly = [];
for i = 1:N
    R_poly_i = 0;
    for j = 1:N
        R_poly_i = R_poly_i + Gam(j) * poly(i, j);
    end
    R_poly = [R_poly, R_poly_i];
end
















%==========================================================================
% forming result
% result.eta = eta;
result.eigenvalues = eig2_LFP';

result.eigenvalues_Rx_ps = eig2';

result.eigenvalues_Rx_sp = eig4';

result.eigenfunction_Rx = eigf1;

result.lambda0_R = eig2(1);
result.lambda0 = eig2_LFP(1);

result.quadpt = pt';
result.weight = wt;
result.Rx = R;
result.Rpoly = R_poly;
result.Gam = Gam;
result.poly = poly;

result.L_ps = L_ps;
result.L_sp = L_sp;
result.D_ps = D_ps;
result.D_sp = D_sp;


% ???
% result.eta = 1 - eig2(1)/ 2 * (1/sqrt(E)) * exp(-E);
result.etahat = eig2(1) / (exp(-E));

% c = linsolve(D, Gam');
% eta = 0;
% for i = 1:N
%     eta = eta - c(i) * Gam(i);
% end
% 
% result.D = D;
% result.c = c;
% result.eta = eta;

end