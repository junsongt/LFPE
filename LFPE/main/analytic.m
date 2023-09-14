function result = analytic(rm, E)
% eps = linspace(0, 10, 1000); 
% eps = [3];
% m_1 = 1; 
% m_2 = 100; 


xm1 = rm/(1+rm);
xm2 = 1 - xm1;

S_00 = 1; 
S_10 = gamma(2.5)/(gamma(1.5)); 
S_11 = -1; 
S_20 = gamma(3.5)/gamma(1.5)/2; 
S_21 = -gamma(3.5)/(gamma(2.5)); 
S_22 = 1/2; 

int_1 = @(x) exp(-x.^2).*x.^3.*(1 - E./x)./4; 
int_2 = @(x) exp(-x.^2).*x.^5.*(1 - E./x)./4; 
int_3 = @(x) exp(-x.^2).*x.^7.*(1 - E./x)./4;
K_0 = 1./pi.*integral( @(x) exp(-x.^2).*x.^3.*(1 - E./x)./4, E, Inf);
K_1 = 1./pi.*integral( @(x) exp(-x.^2).*x.^5.*(1 - E./x)./4, E, Inf);
K_2 = 1./pi.*integral( @(x) exp(-x.^2).*x.^7.*(1 - E./x)./4, E, Inf);


%A_01 = S_00*K_0; 
%A_02 = S_00*K_0; 


%A_11 = (1 - m_1/(m_1 + m_2))*(S_10*K_0 + S_11*K_1); 
%A_21 = (1 - m_2/(m_1 + m_2))*(S_10*K_0 + S_11*K_1); 
%A_12 = (1 - m_1/(m_1 + m_2))^2*(S_20*K_0 + S_21*K_1 + S_22*K_2);

%alpha_12 = A_12; 
%alpha_11 = A_11; 

m_ratio = 1/2; 
%m_ratio = m_1/(m_1 + m_2);

M_11 = -8*m_ratio*(1-m_ratio); 
M_12 = 4*m_ratio*(1-m_ratio)^2; 
M_22 = -2*m_ratio*(1-m_ratio)*(15*m_ratio^2 - 18*m_ratio + 13); 

%a_11 = (M_12*alpha_12 - M_22*alpha_11)/(M_12^2 - M_11*M_22); 
%a_12 = (M_12*alpha_11 - M_11*alpha_12)/(M_12^2 - M_11*M_22); 

%eta = -(a_11*A_11/A_01 + a_12*A_12/A_02)

A_01 = 4 * S_00*K_0;

A_11 = 4 *(1 - rm/(1+rm))*(S_10*K_0 + S_11*K_1);

A_21 = 4 * (1 - rm/(1+rm))^2*(S_20*K_0 + S_21*K_1 + S_22*K_2);

alpha(1) = A_11;
alpha(2) = A_21;


M(1,1) = M_11;
M(1,2) = M_12;
M(2,1) = M_12;
M(2,2) = M_22;


a = linsolve(M, alpha');
etahat = -(a(1)*A_11 + a(2)*A_21) / A_01;




a_11 = (M_12*A_21 - M_22*A_11)/(M_12^2 - M_11*M_22);
a_12 = (M_12*A_11 - M_11*A_21)/(M_12^2 - M_11*M_22);
eta = -(a_11*A_11/A_01 + a_12*A_21/A_01);


result.etahat = etahat;
result.eta = eta;
result.a = a;
result.a11 = a_11;
result.a12 = a_12;
result.M = M;
result.alpha = alpha;





% plot(eps,eta,'LineWidth', 2, 'color', 'black')
% set(gca,'FontSize',24)
% xlabel('$\epsilon/k_bT$','Interpreter','latex','fontsize',24)
% ylabel('$\eta$','Interpreter','Latex','fontsize',24)


end