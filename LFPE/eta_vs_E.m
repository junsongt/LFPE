close; clear; clc;

labels = ["a", "b", "c", "d","e","f"];

% mass_ratios = [1, 1/10, 1/100, 1/1000, 1/10000];
mass_ratios = [0.1, 0.05, 0.02, 0.01];
% mass_ratios = [5,2,1,0.8,0.5];
% mass_ratios = [5, 2, 1];
% mass_ratios = [0.6, 0.3, 0.15, 0.05, 0.02, 0.01];
% mass_ratios = [1, 1.5, 2, 2.5, 3, 4];
% mass_ratios = [1, 0.5, 0.2, 0.1, 0.08, 0.05, 0.01];
% mass_ratios = [0.5, 0.2, 0.1, 0.05];
% mass_ratios = [1, 0.5, 0.1, 0.05, 0.01, 0.005, 0.001];

% n_ratios = [1, 1/10, 1/100, 1/1000, 1/10000];
% n_ratios = [1, 1/2, 1/4, 1/8];
% n_ratios = [1, 1.5, 2, 2.5, 3, 4, 5];
% Es = linspace(0,10,100);

% new E domain: [2,10] for Lorentz limit
Es = linspace(2,10,100);

rd = 1;
rm = 1;
rn = 1;
N = 50;
% Ns = [1,2,3,4];
% N = 65 is bound


% % PAPER II:
% for i = 1 : length(n_ratios)
%     etas = [];
%     rm = 1;
%     rn = n_ratios(i);
%     for j = 1 : length(Es)
%         E = Es(j);
%         result = NEQ2(rm,rn,rd,E,N);
%         eta = result.eta;
%         etas = [etas eta];
%     end
%     plot(Es, etas);
%     title("Fig1, paper II")
%     % title("mass ratio to 0")
%     xlabel('$\varepsilon^* / KT$','Interpreter','latex','fontsize',20);
%     ylabel('$\eta$','Interpreter','Latex','fontsize',20);
%     legend('$m_1/m_2=1$', ...
%         '$m_1/m_2=1.5$', ...
%         '$m_1/m_2=2$', ...
%         '$m_1/m_2=2.5$', ...
%         '$m_1/m_2=3$', ...
%         '$m_1/m_2=4$', 'Interpreter',"latex");
%     % legend('$m_1/m_2=1$', ...
%     %     '$m_1/m_2=1/10$', ...
%     %     '$m_1/m_2=1/100$', ...
%     %     '$m_1/m_2=1/1000$', ...
%     %     'Interpreter',"latex");
%     hold on;
%
% end
% savefig("paper II")
% savefig("massratio")



% % check Sonine convergence of eta
% Ns = [1,2,3,4];
%
% for i = 1:length(Ns)
%     etas = [];
%     rn = 1/10000;
%     n = Ns(i);
%     for j = 1 : length(Es)
%         E = Es(j);
%         result = NEQ2(rm,rn,rd,E,n);
%         eta = result.eta;
%         etas = [etas eta];
%     end
%     plot(Es, etas);
%     title("Different $N$, $\frac{n_1}{n_2} = \frac{1}{10000}, " + ...
%         "\frac{m_1}{m_2} = 1$",'Interpreter','latex')
%     xlabel('$\varepsilon^* / KT$','Interpreter','latex','fontsize',20);
%     ylabel('$\eta$','Interpreter','Latex','fontsize',20);
%     legend('N = 1', 'N = 2', 'N = 3', 'N = 4');
%     hold on;
% end
% savefig("Nconvergence")


% PAPER III:

% % plots of eta vs E for different mass ratios
% for i = 1 : length(mass_ratios)
%     etas = [];
%     rm = mass_ratios(i);
%     rn = 1;
%     delta = 0.01;
%     max_index = 1;
%     max_eta = 0;
%     for j = 1 : length(Es)
%         E = Es(j);
%         result = NEQ3(rm,E,rn,delta, N);
%         eta = result.eta;
%         etas = [etas eta];
%         if eta > max_eta
%             max_eta = eta;
%             max_index = j;
%         end
%     end
% 
%     plot(Es, etas, "black");
%     title("N=10, Sonine basis, dR/dE = 0.01, LOC");
%     str = strcat("(", labels(i), ")");
%     text(Es(max_index),etas(max_index), str,'fontsize',15);
%     % x_start = 0.5;
%     % x_end = 0.25;
%     % y_start = etas(max_index)/0.25;
%     % y_end = etas(max_index)/0.25;
%     % annotation('textarrow',[x_start x_end],[y_start y_end],'String', labels(i), 'fontsize',15);
%     xlabel('$\varepsilon^* / KT$','Interpreter','latex','fontsize',15);
%     ylabel('$\eta$','Interpreter','Latex','fontsize',15);
%     legend('$(a): \frac{m}{M}=0.1$', ...
%         '$(b): \frac{m}{M}=0.05$', ...
%         '$(c): \frac{m}{M}=0.02$', ...
%         '$(d): \frac{m}{M}=0.01$', ...
%         'Interpreter','Latex','fontsize',15);
%     hold on;
% 
% end
% % savefig("paperIII")

% % Convergence of plots of eta vs E for certain choice of mass ratio
% % rn = n_ratios(i);
% rm = 0.1;
% rn = 1;
% delta = 0.01;
% for i = 1 : length(Ns)
%     N = Ns(i);
%     etas = [];
%     max_index = 1;
%     max_eta = 0;
%     for j = 1 : length(Es)
%         E = Es(j);
%         result = NEQ3(rm,E,rn,delta, N);
%         eta = result.eta;
%         etas = [etas eta];
%         if eta > max_eta
%             max_eta = eta;
%             max_index = j;
%         end
%     end
% 
%     plot(Es, etas, "black");
%     title("m/M = 0.1, Sonine basis, LOC");
%     str = strcat("(", labels(i), ")");
%     % text(Es(max_index),etas(max_index), str,'fontsize',20);
%     text(Es(end-10*i),etas(end-10*i), str,'fontsize',20);
%     % x_start = 0.5;
%     % x_end = 0.25;
%     % y_start = etas(max_index)/0.25;
%     % y_end = etas(max_index)/0.25;
%     % annotation('textarrow',[x_start x_end],[y_start y_end],'String', labels(i), 'fontsize',15);
%     xlabel('$\varepsilon^* / KT$','Interpreter','latex','fontsize',20);
%     ylabel('$\eta$','Interpreter','Latex','fontsize',20);
%     legend('$(a): N=1$', ...
%         '$(b): N=2$', ...
%         '$(c): N=3$', ...
%         '$(d): N=4$', ...
%         '$(e): N=6$', ...
%         '$(f): N=8$', ...
%         '$(g): N=10$', ...
%         'Interpreter','Latex','fontsize',15);
%     hold on;
% 
% end








% % convergence of eta in paper II to eta in paper III when n1/n2 -> 0
% etasIII = [];
% rm = 1;
% delta = 1;
% for j = 1 : length(Es)
%     E = Es(j);
%     result = NEQ3(rm,E,delta,N);
%     % result = NEQ4(rm,rn,rd,E,N);
%     eta = result.eta;
%     etasIII = [etasIII eta];
% end
% plot(Es, etasIII);
% title("convergence of $\eta$ to isothermal case when $\frac{n_1}{n_2} \to 0$", 'Interpreter', "latex")
% xlabel('$\varepsilon^* / KT$','Interpreter','latex','fontsize',20);
% ylabel('$\eta$','Interpreter','Latex','fontsize',20);
% hold on;
%
%
% for i = 1:length(n_ratios)
%     etasII = [];
%     rn = n_ratios(i);
%     for j = 1 : length(Es)
%         E = Es(j);
%         result = NEQ2(rm,rn,rd,E,N);
%         eta = result.eta;
%         etasII = [etasII eta];
%     end
%
%     plot(Es, etasII);
%     xlabel('$\varepsilon^* / KT$','Interpreter','latex','fontsize',20);
%     ylabel('$\eta$','Interpreter','Latex','fontsize',20);
%     legend('isothermal', '$\frac{n_1}{n_2}=1$', ...
%         '$\frac{n_1}{n_2}=\frac{1}{10}$', ...
%         '$\frac{n_1}{n_2}=\frac{1}{100}$', ...
%         '$\frac{n_1}{n_2}=\frac{1}{1000}$', ...
%         '$\frac{n_1}{n_2}=\frac{1}{10000}$','Interpreter','latex');
%     hold on;
% end
% savefig("convergenceisothermal")



% Spectral method
for i = 1 : length(mass_ratios)
    etas = [];
    rm = mass_ratios(i);
    rn = 1;
    max_index = 1;
    max_eta = 0;
    for j = 1 : length(Es)
        E = Es(j);
        result = LFPE_Rx(N, 1, E, rm);
        etahat = result.etahat;
        % k = result.lambda1_R;
        % keq = result.lambda1;
        % eta = 1 - k/keq;
        etas = [etas etahat];
        if etahat > max_eta
            max_eta = etahat;
            max_index = j;
        end
    end
    plot(Es, etas);
    % title("N = 40, Maxwell polynomial basis");
    str = strcat("(", labels(i), ")");
    % text(Es(10),etas(10), labels(i));
    % text(Es(end-10*i),etas(end-10*i), str,'fontsize',15);
    text(Es(max_index),etas(max_index), str,'fontsize',15);
    xlabel('$\varepsilon^* / KT$','Interpreter','latex','fontsize',15);
    ylabel('$1 - \frac{\lambda_0}{k_{eq}}$','Interpreter','Latex','fontsize',15);
    legend('$(a): \frac{m}{M}=0.1$', ...
        '$(b): \frac{m}{M}=0.05$', ...
        '$(c): \frac{m}{M}=0.02$', ...
        '$(d): \frac{m}{M}=0.01$', ...
        'Interpreter','Latex','fontsize',10);
    hold on;

end
savefig("etahat")




