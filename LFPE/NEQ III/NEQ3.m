% Author: Junsong Tang
% This function is based on Shizgal, 1971(Paper III), which produces an object containing:
% eta
% Main matrix: AA collision matrix
% alpha quantity
% A integral(moment)
% coefficient vector: a


function res = NEQ3(rm, E, rn, delta, N)
%RM is the mass ratio: n/M
%RN is the number density ratio: n/n_M
%E = E*/kt the reduced threshold energy in the line of centers cross
%section
%N is the number of Sonine Polynomial expansion retained, note the program
%will not work if N is greater than 10.
% delta is (dR/dE)^2


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
xm2 = 1 - xm1;

% AA collision  matrix
a11 = mat1(xm1);
% a11 = mat1(1/2);


% Set up the main matrix
% Divided eqs. by Q12*N1*N2; Factor in front of A(I,J)
% is N1*N1*Q/(N1*N2*Q12)

fac1 = 1 * sqrt(xm2 / 2) * (2 * (1/2))^2;
%fac1 = (xn1 / xn2) * sqrt(xm2 / 2) * (2 * d1 / (d1 + d2))^2;



% Set up the inhomogeneous column vector
result = alf3(E, N, xm2);
alpha = result.alpha;
A = result.A;

for i = 1:N
   rhs(i) = alpha(i);
end



% Solve the set of simultaneous equations
mtx = a11(1:N, 1:N);

% some scalar: (n1/n2) * (dR/dE)^2 before alpha quantity
rhs = rhs * rn * delta;

a = linsolve(mtx, rhs');
% a = rn * a;

sum = 0;
for i = 1:N
    % sum = sum - rhs(i) * xm2^i * aa(i+1) / aa(1); % index!!!
    sum = sum - a(i) * A(i+1) / A(1);
end





% The below is based on Shizgal, 1974(time dependent isothermal)
% Now use (1 - lambda_0/keq) to estimate eta

% elastic matrix
Elastic(2:N+1, 2:N+1) = mtx;
for i = 1:N+1
    Elastic(1,i) = 0;
    Elastic(i,1) = 0;
end

% reaction matrix
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



% main matrix M = E - R ===> normalizing M
for i = 1:N+1
    for j = 1:N+1
        Ni = 2 * (gamma(i-1+1.5)/factorial(i-1)) / sqrt(pi);
        Nj = 2 * (gamma(j-1+1.5)/factorial(j-1)) / sqrt(pi);
        fac = 1 / (sqrt(Ni) * sqrt(Nj));
        M(i,j) = fac * (Elastic(i,j) - R(i,j));
    end
end

% calculate eigenvalues of normalized M
[eigf,eigv]=eig(M);
eigen1 = diag(eigv);
eigen2 = sort(abs(eigen1),'ascend')';


% Forming results
res.eta = sum;
res.mtx = mtx;
res.rhs = rhs;
res.a = a';
res.A = A;
res.alpha = alpha;

res.eigen_values = eigen2;
res.Elastic = Elastic;
res.Reaction = R;
res.M = M;
res.etahat = 1 - abs(eigen2(1))/A(1);
end



%==========================================================================
% Helper functions

function result = alf3(E, n, xm2)
% xm2 := m2/(m1+m2)
% WARNING: 
% s has 0 based index;
% xk has 0 based index;
% xj has 0 based index;
% a has 0 based index;
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
end

% Calculate the K(k) integrals

% Uncomment below if for LOC cross section case:
xk(1) = expe / 2;
for k = 2:nn+1
    xk(k) = xj(k+1) - E * xj(k);
end


% % Uncomment below if for STEP cross section case:
% xk(1)=(1/2) *(E + 1) * expe;
% for k=2:nn+1
%     xk(k) = 0.5*E^k*expe + k*xk(k-1);
% end


% %  Calculate K(k) integral (explicit integral solving without using recursive relation)
% for k = 1:nn+1
%     int_k = @(x) exp(-x.^2).*x.^(2*(k-1)+3).*(1 - E./x)./4;
%     xk(k) = 1./pi.*integral( @(x) exp(-x.^2).*x.^(2*(k-1)+3).*(1 - E./x)./4, E, Inf);
% end



% Calculate the A(j) integrals
% a = zeros(nn+1, 1);
for i = 1:nn+1
    sum = 0;
    for k = 1:i+1
        sum = sum + s(i, k) * xk(k);
    end
    % a(i) = 4 * sum;
    A(i) = 4 * (xm2^(i-1)) * sum;
end


% Calculate the Alpha(i) quantities
% alpha(1) = xm2 * a(2); % index!!!
alpha(1) = A(2); % index!!!
for k = 2:nn
    % alpha(k) = xm2^k * a(k+1); % index!!!
    alpha(k) = A(k+1);
end

result.alpha = alpha;
result.A = A;
result.s = s;
result.K = xk;
end




% Table II matrix extended to 10x10
function a11 = mat1(xm)
ym=-1+xm;
xm2=xm*xm;
xm3=xm*xm2;
xm4=xm*xm3;
xm5=xm*xm4;
xm6=xm*xm5;
xm7=xm*xm6;
xm8=xm*xm7;
xm9=xm*xm8;
xm10=xm*xm9;
xm11=xm*xm10;
xm12=xm*xm11;
xm13=xm*xm12;
xm14=xm*xm13;
xm15=xm*xm14;
xm16=xm*xm15;
xm17=xm*xm16;
xm18=xm*xm17;

a11(1,1)=8*xm*ym ;
a11(1,2)=4*xm*ym^2;
a11(1,3)=-xm*ym^3;

% !!!
a11(1,4)=xm*ym^4/2;  % from code
% a11(1,4)=-xm*ym^4/2; % from paper

a11(1,5)=-5*xm*ym^5/16;
a11(1,6)=7*xm*ym^6/32;
a11(1,7)=-21*xm*ym^7/128;
a11(1,8)=33*xm*ym^8/256;
a11(1,9)=-429*xm*ym^9/4096;
a11(1,10)=715*xm*ym^10/8192;
a11(2,2)=2*xm*ym*(15*xm2-18*xm+ 13);
a11(2,3)=xm*(35*xm2-30*xm+ 23)*ym^2/2;
a11(2,4)=-xm*(21*xm2-14*xm+ 11)*ym^3/4;
a11(2,5)=xm*(99*xm2-54*xm+ 43)*ym^4/32;
a11(2,6)=-xm*(143*xm2-66*xm+ 53)*ym^5/64;
a11(2,7)=7*xm*(65*xm2-26*xm+ 21)*ym^6/256;
a11(2,8)=-3*xm*(255*xm2-90*xm+ 73)*ym^7/512;
a11(2,9)=33*xm*(323*xm2-102*xm+ 83)*ym^8/8192;
a11(2,10)=-143*xm*(133*xm2-38*xm+ 31)*ym^9/16384;

a11(3,3)=xm*ym*( 945*xm4-2100*xm3+2190*xm2-1188*xm+ 433)/8;
a11(3,4)=xm*(1155*xm4-2100*xm3+2030*xm2-940*xm+ 359)*ym^2/16; % from code
% a11(3,4)=-xm*(1155*xm4-2100*xm3+2030*xm2-940*xm+ 359)*ym^2/16; % from paper
a11(3,5)=-xm*(3003*xm4-4620*xm3+4242*xm2-1708*xm+ 667)*ym^3/128;
a11(3,6)=xm*(3861*xm4-5148*xm3+4554*xm2-1620*xm+ 641)*ym^4/256;
a11(3,7)=-xm*(12155*xm4-14300*xm3+12298*xm2-3916*xm+ 1563)*ym^5/1024;
a11(3,8)=xm*(20995*xm4-22100*xm3+18590*xm2-5356*xm+2151)*ym^6/2048;
a11(3,9)=-3*xm*(101745*xm4-96900*xm3+80070*xm2-21060*xm+ 8497)*ym^7/32768;
a11(3,10)=11*xm*(52003*xm4-45220*xm3+36822*xm2-8908*xm+3607)*ym^8/65536;

a11(4,4)=xm*ym*( 15015*xm6-48510*xm5+72135*xm4-62300*xm3 +34485*xm2-12102*xm+ 2957)/32;
a11(4,5)=xm*(75075*xm6-210210*xm5+285285*xm4-224700*xm3+117005*xm2-37090*xm+ 9419)*ym^2/256;
a11(4,6)=-xm*(51051*xm6-126126*xm5+159159*xm4-115500*xm3+57393*xm2-16534*xm+ 4285)*ym^3/512;
a11(4,7)=xm*(138567*xm6-306306*xm5+364221*xm4-245388*xm3+117513*xm2-30978*xm+ 8131)*ym^4/2048;
a11(4,8)=-xm*(230945*xm6-461890*xm5+522665*xm4-328900*xm3 +152867*xm2-37114*xm+ 9827)*ym^5/4096;
a11(4,9)=xm*(3380195*xm6-6172530*xm5+6697405*xm4-3955900*xm3+1793805*xm2-403442*xm+ 107507)*ym^6/65536;
a11(4,10)=- xm*(6500375*xm6-10920630*xm5+11429355*xm4-6363100*xm3+2826165*xm2-591870*xm+ 158489)*ym^7/131072;

a11(5,5)=xm*ym*( 3828825*xm8-16216200*xm7+31771740*xm6-37449720*xm5+29343510*xm4-15894200*xm3+6043740*xm2-1568328*xm+288473)/2048;
a11(5,6)=xm*(4849845*xm8-18378360*xm7+33093060*xm6-36156120*xm5+26539590*xm4-13491240*xm3+4879700*xm2-1175560*xm+223469) *ym^2/4096;
a11(5,7)=-xm*(6789783*xm8-23279256*xm7+39002964*xm6-39855816*xm5 +27669642*xm4-13281576*xm3+4611348*xm2-1033592*xm+200183)*ym^3 /16384;
a11(5,8)=xm*(9561123*xm8-29930472*xm7+47112780*xm6-45333288*xm5+29989674*xm4-13655928*xm3+4582908*xm2-958968*xm+ 188011)* ym^4/32768;
a11(5,9)=-xm*(132793375*xm8-382444920*xm7+569972260*xm6 -519164360*xm5+329230330*xm4-142742600*xm3+46551076*xm2- 9126392*xm+ 1804831)*ym^5/524288;
a11(5,10)=xm*(253514625*xm8-676039000*xm7+959975380*xm6- 831234040*xm5+507785070*xm4-210259400*xm3+66911780*xm2 -12333672*xm+2454937)*ym^6/1048576;

a11(6,6)=xm*ym*(61108047*xm10-320089770*xm9+781846065*xm8 -1170809640*xm7 +1194383190*xm6-872733708*xm5+467955810*xm4-185239880*xm3 +53790315*xm2 -11075802*xm+1634141)/8192;
a11(6,7)=xm*(156165009*xm10-746876130*xm9+1692595905*xm8 -2370808440*xm7 +2277445170*xm6-1573968396*xm5+801955770*xm4-302163960*xm3+84127645*xm2 -16341410*xm+ 2481445)*ym^2/32768;
a11(6,8)=-xm*(111546435*xm10-490804314*xm9+1040776737*xm8 -1373476104*xm7 +1251055806*xm6-822593772*xm5+400642242*xm4-144297384*xm3+38760687*xm2-7107338*xm+ 1098029)*ym^3/65536;
a11(6,9)=xm*(1290751605*xm10-5258617650*xm9+10507674177*xm8 -13139477208*xm7 +11412100986*xm6-7171848684*xm5+3355288794*xm4-1158996696*xm3 +301884561*xm2 -52347186*xm+ 8182373)*ym^4/1048576;
a11(6,10)=-xm*(2310604725*xm10-8764362750*xm9+16599171875*xm8 -19759654200*xm7 +16440050770*xm6-9911235620*xm5+4472091910*xm4-1485512600*xm3 +376733929*xm2 -61909518*xm+ 9759719)*ym^5/2097152;

a11(7,7)=xm*ym*(3904125225*xm12-24361741404*xm11+71414937594*xm10 -129956446620*xm9+163735106535*xm8-150911200440*xm7+104843054316*xm6 -55735157016*xm5+22772818215*xm4-7112122220*xm3+1672265850*xm2 -285428940*xm+35164265)/131072;
a11(7,8)=xm*(94408854540*xm6+162649251765*xm8-136358242020*xm9 -29002073100*xm11 +1278516050*xm2-5638029180*xm3+18802525665*xm4- 47993873928*xm5+26292995 +5019589575*xm12-208028500*xm+79554917442*xm10 -142487425080*xm7)*ym^2/262144;
a11(7,9)=-xm*(396053453796*xm6+745573461633*xm8-656696172132*xm9-156611194740*xm11+4626970110*xm2-21061507236*xm3+72925101249*xm4 -193337248104*xm5+29113619535*xm12-717680180*xm+ 404467373310*xm10 -623860781544*xm7+92174735)*ym^3/4194304;
a11(7,10)=xm*(434414750484*xm6+886631619159*xm8-816838608300*xm9-216272602260*xm11+4436214750*xm2-20764720548*xm3+74457949995*xm4-204323430168*xm5+42977247885*xm12-656558460*xm+ 528634490670*xm10 -711457296264*xm7+85274065)*ym^4/8388608;


%9,10 element
f=42854309105;
a910(1)=-452887954000/f;
a910(2)=3687635894600/f;
a910(3)=-22245723005520/f;
a910(4)=104265844271820/f;
a910(5)=-387650968937232/f;
a910(6)=1159665950603160/f;
a910(7)=-2812813034837520/f;
a910(8)=5546913045869910/f;
a910(9)=-8875607053513200/f;
a910(10)=11443235187311928/f;
a910(11)=-11733194701628400/f;
a910(12)=9363629399424300/f;
a910(13)=-5616987662286000/f;
a910(14)=2387815892490600/f;
a910(15)=-642939628359600/f;
a910(16)=82731202178625/f;
sum1=0.0;
for j=1:16
    sum1=sum1+a910(j)*xm^j;
end
a11(9,10)=f*xm*ym^2*(1+sum1)/268435456;

%8,10 element
f=467774345;
a810(1)=-4257776390/f;
a810(2)=31684909305/f;
a810(3)=-170243480484/f;
a810(4)=705459111357/f;
a810(5)=-2286626728386/f;
a810(6)=5872296418989/f;
a810(7)=-11990037338136/f;
a810(8)=19419597063867/f;
a810(9)=-24709543188330/f;
a810(10)=24246737029515/f;
a810(11)=-17766669981060/f;
a810(12)=9190199233215/f;
a810(13)=-3008407351950/f;
a810(14)=472749726735/f;
sum2=0.0;
for j=1:14
    sum2=sum2+a810(j)*xm^j;
end
a11(8,10)=-f*xm*ym^3*(1+sum2)/16777216;


% 8,9 element
f=1076028115;
a89(1)=-9942043850/f;
a89(2)=71026222225/f;
a89(3)=-370954581060/f;
a89(4)=1488337907595/f;
a89(5)=-4668724706646/f;
a89(6)=11583779512305/f;
a89(7)=-22815544597560/f;
a89(8)=35577827590305/f;
a89(9)=-43490917139670/f;
a89(10)=40893748514619/f;
a89(11)=-28625046149700/f;
a89(12)=14089987937025/f;
a89(13)=-4367042930250/f;
a89(14)=644658718275/f;
sum3=0.0;
for j=1:14
    sum3=sum3+a89(j)*xm^j;
end
a11(8,9)=f*xm*ym^2*(1+sum3)/8388608;


%8,8 element
f=182066285;
a8(1)=-1721444910/f;
a8(2)=11810004225/f;
a8(3)=-59753862980/f;
a8(4)=231559645905/f;
a8(5)=-700830690714/f;
a8(6)=1674775691589/f;
a8(7)=-3171348025560/f;
a8(8)=4744382057415/f;
a8(9)=-5550676701570/f;
a8(10)=4981158433251/f;
a8(11)=-3316677079716/f;
a8(12)=1546591321275/f;
a8(13)=-451763061750/f;
a8(14)=62386327575/f;
sum4=0.0;
for j=1:14
    sum4=sum4+a8(j)*xm^j;
end
a11(8,8)=f*xm*ym*(1+sum4)/524288;

%9,9 element
f=58576234985;
a9(1)=-632166153360/f;
a9(2)=4970036355000/f;
a9(3)=-29147104802480/f;
a9(4)=132528511237500/f;
a9(5)=-477592248372624/f;
a9(6)=1383104916737544/f;
a9(7)=-3243422074566960/f;
a9(8)=6174719656638390/f;
a9(9)=-9522746381594160/f;
a9(10)=11811837017916936/f;
a9(11)=-11627608594217616/f;
a9(12)=8887700917995900/f;
a9(13)=-5092273232046000/f;
a9(14)=2060745172457400/f;
a9(15)=-526041514112400/f;
a9(16)=63821213109225/f;
sum5=0.0;
for j=1:16
    sum5=sum5+a9(j)*xm^j;
end
a11(9,9)=f*xm*ym*(1+sum5)/134217728;


% 10, 10 element
f=287640532965;
a10(1)=-3488705287890/f;
a10(2)=30920613176175/f;
a10(3)=-206192215491920/f;
a10(4)=1075801164307380/f;
a10(5)=-4496334554189016/f;
a10(6)=15291169120022796/f;
a10(7)=-42734949206217840/f;
a10(8)=98697853772485110/f;
a10(9)=-188727042254679500/f;
a10(10)=298311442465070034/f;
a10(11)=-387725926434550704/f;
a10(12)=410419720867256100/f;
a10(13)=-348460912214415000/f;
a10(14)=231809210097507900/f;
a10(15)=-116488970847334800/f;
a10(16)=41618522193115725/f;
a10(17)=-9431357048363250/f;
a10(18)=1020351493536375/f;
sum6=0.0;
for j=1:18 
    sum6=sum6+a10(j)*xm^j;
end
a11(10,10)=f*xm*ym*(1+sum6)/536870912;


for i = 1:9
    for j = i+1:10
        a11(j, i) = a11(i, j);
    end
end

end
