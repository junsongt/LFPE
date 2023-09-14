

eps = linspace(0, 10, 1000); 
% eps = [3];
m_1 = 1; 
m_2 = 100; 

S_00 = 1; 
S_10 = gamma(2.5)./(gamma(1.5)); 
S_11 = -1; 
S_20 = gamma(3.5)./gamma(1.5)./2; 
S_21 = -gamma(3.5)./(gamma(2.5)); 
S_22 = 1./2; 

int_1 = @(x) exp(-x.^2).*x.^3.*(1 - eps./x)./4; 
int_2 = @(x) exp(-x.^2).*x.^5.*(1 - eps./x)./4; 
int_3 = @(x) exp(-x.^2).*x.^7.*(1 - eps./x)./4;

for i= 1:length(eps)
    K_0(i) = 1./pi.*integral( @(x) exp(-x.^2).*x.^3.*(1 - eps(i)./x)./4, eps(i), Inf); 
    K_1(i) = 1./pi.*integral( @(x) exp(-x.^2).*x.^5.*(1 - eps(i)./x)./4, eps(i), Inf); 
    K_2(i) = 1./pi.*integral( @(x) exp(-x.^2).*x.^7.*(1 - eps(i)./x)./4, eps(i), Inf); 
end

%A_01 = S_00.*K_0; 
%A_02 = S_00.*K_0; 


%A_11 = (1 - m_1./(m_1 + m_2)).*(S_10.*K_0 + S_11.*K_1); 
%A_21 = (1 - m_2./(m_1 + m_2)).*(S_10.*K_0 + S_11.*K_1); 
%A_12 = (1 - m_1./(m_1 + m_2)).^2.*(S_20.*K_0 + S_21.*K_1 + S_22.*K_2);

%alpha_12 = A_12; 
%alpha_11 = A_11; 

m_ratio = 1./2; 
%m_ratio = m_1./(m_1 + m_2);

M_11 = -8.*m_ratio.*(1-m_ratio); 
M_12 = 4.*m_ratio.*(1-m_ratio).^2; 
M_22 = -2.*m_ratio.*(1-m_ratio).*(15.*m_ratio.^2 - 18.*m_ratio + 13); 

%a_11 = (M_12.*alpha_12 - M_22.*alpha_11)./(M_12.^2 - M_11.*M_22); 
%a_12 = (M_12.*alpha_11 - M_11.*alpha_12)./(M_12.^2 - M_11.*M_22); 

%eta = -(a_11.*A_11./A_01 + a_12.*A_12./A_02)

for i=1:length(eps)
    A_01(i) = 4*S_00.*K_0(i); 
    A_02(i) = 4*S_00.*K_0(i); 
    A_11(i) = 4*(1 - m_1./(m_1 + m_2)).*(S_10.*K_0(i) + S_11.*K_1(i)); 
    A_21(i) = 4*(1 - m_2./(m_1 + m_2)).*(S_10.*K_0(i) + S_11.*K_1(i)); 
    A_12(i) = 4*(1 - m_1./(m_1 + m_2)).^2.*(S_20.*K_0(i) + S_21.*K_1(i) + S_22.*K_2(i));
end

for i=1:length(eps)
    a_11(i) = (M_12.*A_12(i) - M_22.*A_11(i))./(M_12.^2 - M_11.*M_22); 
    a_12(i) = (M_12.*A_11(i) - M_11.*A_12(i))./(M_12.^2 - M_11.*M_22);
    eta(i) = -(a_11(i).*A_11(i)./A_01(i) + a_12(i).*A_12(i)./A_01(i));
end


plot(eps,eta,'LineWidth', 2, 'color', 'black')
set(gca,'FontSize',24)
xlabel('$\epsilon/k_bT$','Interpreter','latex','fontsize',24)
ylabel('$\eta$','Interpreter','Latex','fontsize',24)








