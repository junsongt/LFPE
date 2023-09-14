function ab=ab_maxwell_p2(m,nint,npts,xmin,xmax)
%Guass Stieltjes procedure for caluation of recurrence coefficients
% alpha and betta for weight function w(x)=x*x*exp(-x*x) 
format long e
%nint=# of intervals; npts=# of pts per interval [xmin,xmax]
pwmd = mdquad(nint,npts,xmin,xmax); ntot=nint*npts;
%multi-domain matrix of quad pts p and wts w: 
%weight function w(x)=x*x*exp(-x*x)
p=pwmd(:,1); w=pwmd(:,2); psq=p.*p; wtfcn=psq.*exp(-psq);
b0=ones(ntot,1);
%First two integrals
s1=sum(w.*wtfcn); s2=sum(w.*(p.*wtfcn));
h(1)=s1; alfa(1)=s2/s1; beta(1)=0;k=1;
%Norm and alpha_1 and betta_1
% fprintf('%i %1.16e %1.16e\n',k,alfa(k),beta(k))
%Polynomial P_1(x)
b1=p-alfa(1);
%Next two integrals 
s1=sum(w.*(wtfcn.*(b1.^2))); s2=sum(w.*(p.*(wtfcn.*(b1.^2))));
alfa(2)=s2/s1; h(2)=s1; beta(2)=h(2)/h(1);k=2;
%alpha_2, norm and betta_2
% fprintf('%i %1.16e %1.16e\n',k,alfa(k),beta(k))
for k=3:m
    pma=p-alfa(k-1);
    %Recurrence for the next polynomial
    b2=pma.*b1-beta(k-1)*b0;
    s1=sum(w.*(wtfcn.*(b2.^2)));
    s2=sum(w.*(p.*(wtfcn.*(b2.^2))));
    alfa(k)=s2/s1; h(k)=s1; beta(k)=h(k)/h(k-1);
    %alfa_kl, norm and betta_k
    % fprintf('%i %1.16e %1.16e\n',k,alfa(k),beta(k))
    b0=b1; b1=b2;
end
ab = [alfa; beta]';
end

function pwmd = mdquad(nint,npts,xmin,xmax)
format long e        
ntot=nint*npts;
%Quadrature for mutlidomain with nint intervals and npts per interval
dx=(xmax-xmin)/ntot;
for i=1:nint
    a=xmin+(i-1)*npts*dx; b=a+npts*dx; pw=fejer2(a,b,npts);
    %Fejer quadrature used for each interval
    if i==1
        pwmd=pw;
    else
        pwmd=cat(1,pwmd,pw);
    end
end
end
% FEJER  Fejer quadrature rule.
%        The call pw=fejer(n) generates the n-point Fejer quadrature rule.
function pw=fejer2(a,b,N)
format long e
n=N:-1:1; m=1:floor(N/2); th=(2*n-1)*pi./(2*N); p=cos(th');
for k=N:-1:1
  s=sum(cos(2*m*th(k))./(4*(m.^2)-1));
  w(k)=2*(1-2*s)/N;
end
r1=(b-a)/2.; r2=(a+b)/2.; ps=r1*p+r2; ws=r1*w;
%Map [a,b] onto [-1,1]
pw=[ps,ws'];
end