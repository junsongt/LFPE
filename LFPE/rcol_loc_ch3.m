function rcol_loc_ch3(eps,g)
%eps=5;
m2=1/(1+g);
m1=1-m2;
rtm1=sqrt(m1);
rtm2=sqrt(m2);
h=0.05;
x=[0.00001:h:20];   
nx=length(x);
n=nx-1;
%Calculate gp
rtx=sqrt(x);
t0=(sqrt(eps)+rtm2*rtx)/rtm1;
i0=sqrt(pi)*erfc(t0)/2;
i1=0.5*exp(-t0.^2);
i2=0.5*(exp(-t0.^2).*t0+i0);
g1=m1*i2;
g2=2*rtm1*rtm2*rtx.*i1;
g3=(m2*x-eps).*i0;
gp=rtm1*(g1-g2+g3);
%Calculate gm
t0=(sqrt(eps)-rtm2*rtx)/rtm1;
i0=sqrt(pi)*erfc(t0)/2;
i1=0.5*exp(-t0.^2);
i2=0.5*(exp(-t0.^2).*t0+i0);
g1=m1*i2;
g2=-2*rtm1*rtm2*rtx.*i1;
g3=(m2*x-eps).*i0;
gm=rtm1*(g1-g2+g3);
rcol=(gm-gp)./(rtm1*rtm2*x);
rcol=exp(-x).*sqrt(x).*rcol;
rate=h/3*(rcol(1)+4*sum(rcol(2:2:n))+2*sum(rcol(3:2:n-1))+rcol(n+1));
test=rate/(exp(-eps));
ratio=rate/exp(-eps);
plot(x,1000*rcol/sqrt(2),'-k','linewidth',1.6)
%gtext('$m_1/m_2 = 10$','Interpreter','Latex','fontsize',20)
axis([0 10 0 6])
set(gca,'FontSize',24)
set(gca,'Ytick',[0:2:6],'linewidth',1.6)
set(gca,'Xtick',[0:2:10],'linewidth',1.6)
%str2=num2str(npoly);
%str1={'$n = $'};
%str=strcat(str1,str2);
%text(3,.3,0,str,'Interpreter','latex','fontsize',26)
xlabel('$x$','Interpreter','latex','fontsize',30)
ylabel('$R(x)$','Interpreter','Latex','fontsize',30)