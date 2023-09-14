function rcol_a_spd2(eps)
format long e

%This was used to plot R(y) versus y currently Fig 3.12
%k(T) calcd with speed quad w(x)=x^2exp(x^2)
%gam12=m1/m2
gam12=[.0001 .05 1 10 100];
%fopen(myfile,'%8.4f & %8.4f & %8.4f & %8.4f & %8.4f\n',gam12(1:1:5))
rcolst=[];
nvec=[2:2:20];
load abspeed.dat
for npts=2:2:20
    accst1=[];
    ptwt=abspeed2(npts);
    x=ptwt(:,1);
    w=ptwt(:,2);
    for kk=1:5
        m2=1/(1+gam12(kk));
        m1=1-m2;
        rtm1=sqrt(m1);
        rtm2=sqrt(m2);
        %Calculate gp
        gp=gox(x,eps,m1);
        gm=gox(-x,eps,m1);
        rcol=(gm-gp)./(rtm1*rtm2*x);
        rcolx=2*x.^2.*exp(-x.^2).*rcol/sqrt(pi);
        rcolst=[rcolst rcolx'];
        s=0;
        %for i=1:npts
        %    s=s+w(i)*rcol(i);
        %    fprintf('%i %13.5e %13.5e %9.2f\n',i,w(i),rcol(i),s)
        %    pause
        %end
        rate=2*sum(w.*rcol)/sqrt(pi);
        %rate=2*s/sqrt(pi);
        plot(x,rcolx);
        hold on
        ratio=rate/exp(-eps);
        acc=log10(abs(1-rate/exp(-eps)));

        %fprintf(myfile,'%i %9.2f\n',npts,acc)
        accst1=[accst1 acc];
        %Store each acc value vs N in a vector
    end
    fprintf(myfile,'%i& %8.3f &%8.3f& %8.3f& %8.3f& %8.3f\n',npts,accst1(1),accst1(2),accst1(3),accst1(4),accst1(5))
    %Store each accst1 vector as a column of this matrix
end

%xsq=x.^2;
%plot(xsq,rcolst(:,1),'-k',xsq,rcolst(:,2),'-k',xsq,rcolst(:,3),'-k',xsq,rcolst(:,4),'-k',xsq,rcolst(:,5),'-k','linewidth',1.6)
%xlabel('${y}$','Interpreter','latex','fontsize',26)
%ylabel('$R(y)$','Interpreter','latex','fontsize',26)
%axis([0 15 0 .014])
%set(gca,'FontSize',26)
%set(gca,'Ytick',[0:.002:.014],'linewidth',1.6)
%set(gca,'Xtick',[0:2.5:15],'linewidth',1.6)
%set(gcf, 'Units','centimeters','Papersize',[36,36])
%set(gcf, 'Units','centimeters','Position',[3 3 24 20])
%str1={'$\longleftarrow \frac{m_1}{m_2} = 10^{-4}$'}
%text(6.9,.011,str1,'Interpreter','latex','fontsize',28)
%str2={'$\downarrow \frac{m_1}{m_2} = 100$'}
%text(0.9,.006,str2,'Interpreter','latex','fontsize',28)
function g=gox(x,eps,m1)
m2=1-m1;rtm1=sqrt(m1);rtm2=sqrt(m2);
t0=(sqrt(eps)+rtm2*x)/rtm1;
i0=sqrt(pi)*erfc(t0)/2;
i1=0.5*exp(-t0.^2);
i2=0.5*(exp(-t0.^2).*t0+i0);
g1=m1*i2;
g2=2*rtm1*rtm2*x.*i1;
g3=(m2*x.^2-eps).*i0;
g=rtm1*(g1-g2+g3);
function ptwt = abspeed2(m)
load abspd.dat; 
n=90;
a=abspd(1:m,1);
b=abspd(1:m,2);
rtb=sqrt(b);
rtb(m)=[];
t=diag(rtb,-1)+diag(a)+diag(rtb,1);
[f,lambda]=eig(t);
pt=diag(lambda);
for i=1:m
    wt(i)=sqrt(pi)*f(1,i)^2/4; 
end
ptwt=[pt,wt'];