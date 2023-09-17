% Author: Junsong Tang
% This script could generate:

% (1) the plot of eta ~ E*/KT with different mass ratios in Paper II case

% (2) the plot of convergence of eta ~ E*/KT in increasing number of Sonine
% basis in Paper II case.

% (3) the plot of eta ~ E*/KT with different mass ratios in Paper III case

% (4) the plot of convergence of eta ~ E*/KT in increasing number of Sonine
% basis in Paper III case.

% (5) the plot of eta ~ E*/KT using spectral method in Maxwell polynomial




close; clear; clc;

labels = ["a", "b", "c", "d", "e", "f", "g"];

Es = linspace(0,10,100);

% % new E domain: [2,10] for Lorentz limit
% Es = linspace(2,10,100);


%==========================================================================
% (1) eta ~ E*/KT PAPER II:
N = 4;
rn = 1;
rd = 1;
mass_ratios = [1, 1.5, 2, 2.5, 3, 4];
for i = 1 : length(mass_ratios)
    etas = [];
    rm = mass_ratios(i);
    for j = 1 : length(Es)
        E = Es(j);
        result = NEQ2(rm,rn,rd,E,N);
        eta = result.eta;
        etas = [etas eta];
    end
    plot(Es, etas, "black");
    title("Fig1, paper II");
    str = strcat("(", labels(i), ")");
    text(Es(30),etas(30), str,'fontsize',10);
    xlabel('$\varepsilon^* / KT$','Interpreter','latex','fontsize',15);
    ylabel('$\eta$','Interpreter','Latex','fontsize',15);
    legend('$(a): m_1/m_2=1$', ...
        '$(b): m_1/m_2=1.5$', ...
        '$(c): m_1/m_2=2$', ...
        '$(d): m_1/m_2=2.5$', ...
        '$(e): m_1/m_2=3$', ...
        '$(f): m_1/m_2=4$', ...
        'Interpreter','latex', 'fontsize', 10);
    hold on;
end
savefig("paper II")


%==========================================================================
% % (2) Sonine convergence of eta in Paper II case
% Ns = [1,2,3,4];
% rm = 1/2;
% rn = 1;
% rd = 1;
% for i = 1:length(Ns)
%     etas = [];
%     n = Ns(i);
%     for j = 1 : length(Es)
%         E = Es(j);
%         result = NEQ2(rm,rn,rd,E,n);
%         eta = result.eta;
%         etas = [etas eta];
%     end
%     plot(Es, etas, "black");
%     str = strcat("(", labels(i), ")");
%     text(Es(end-10*i),etas(end-10*i), str,'fontsize',10);
%     xlabel('$\varepsilon^* / KT$','Interpreter','latex','fontsize',15);
%     ylabel('$\eta$','Interpreter','Latex','fontsize',15);
%     legend('$(a): N=1$', ...
%         '$(b): N=2$', ...
%         '$(c): N=3$', ...
%         '$(d): N=4$', ...
%         'Interpreter','Latex','fontsize',10);
%     hold on;
% end
% savefig("Sonine convergence II")




%==========================================================================
% % (3) eta ~ E*/KT PAPER III:
% N = 10;
% mass_ratios = [0.1, 0.05, 0.02, 0.01];
% rn = 1;
% delta = 0.01;
% 
% for i = 1 : length(mass_ratios)
%     etas = [];
%     rm = mass_ratios(i);
% 
%     max_index = 1;
%     max_eta = 0;
%     for j = 1 : length(Es)
%         E = Es(j);
%         result = NEQ3(rm,E,rn,delta,N);
%         eta = result.eta;
%         etas = [etas eta];
% 
%         % update the maximum of eta and the index where eta hits max
%         if eta > max_eta
%             max_eta = eta;
%             max_index = j;
%         end
%     end
% 
%     plot(Es, etas, "black");
%     title("N=10, Sonine basis, (dR/dE)^2 = 0.01, LOC");
%     str = strcat("(", labels(i), ")");
% 
%     % attach labels at each eta~E curve maximum position
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
%         'Interpreter','Latex','fontsize',10);
%     hold on;
% end
% savefig("paper III")


%==========================================================================
% % (4) Sonine convergence of eta in Paper III case
% Ns = [1,2,3,4];
% rm = 0.1;
% rn = 1;
% delta = 0.01;
% for i = 1 : length(Ns)
%     N = Ns(i);
%     etas = [];
%     for j = 1 : length(Es)
%         E = Es(j);
%         result = NEQ3(rm,E,rn,delta, N);
%         eta = result.eta;
%         etas = [etas eta];
%     end
% 
%     plot(Es, etas, "black");
%     title("m/M = 0.1, (dR/dE)^2 = 0.01, Sonine basis, LOC");
%     str = strcat("(", labels(i), ")");
%     % attach labels at the tail of each curve
%     text(Es(end-10*i),etas(end-10*i), str,'fontsize',15);
%     % x_start = 0.5;
%     % x_end = 0.25;
%     % y_start = etas(max_index)/0.25;
%     % y_end = etas(max_index)/0.25;
%     % annotation('textarrow',[x_start x_end],[y_start y_end],'String', labels(i), 'fontsize',15);
%     xlabel('$\varepsilon^* / KT$','Interpreter','latex','fontsize',15);
%     ylabel('$\eta$','Interpreter','Latex','fontsize',15);
%     legend('$(a): N=1$', ...
%         '$(b): N=2$', ...
%         '$(c): N=3$', ...
%         '$(d): N=4$', ...
%         'Interpreter','Latex','fontsize',10);
%     hold on;
% end
% savefig("Sonine convergence III")







%==========================================================================
% % (5) eta ~ E*/KT by pseudospectral/spectral method
% N = 50;
% mass_ratios = [0.1, 0.05, 0.02, 0.01];
% for i = 1 : length(mass_ratios)
%     etas = [];
%     rm = mass_ratios(i);
%     max_index = 1;
%     max_eta = 0;
%     for j = 1 : length(Es)
%         E = Es(j);
%         result = LFPE_Rx(N, 1, E, rm);
%         etahat = result.etahat;
%         etas = [etas etahat];
%         % update the maximum of eta and the index where eta hits max
%         if etahat > max_eta
%             max_eta = etahat;
%             max_index = j;
%         end
%     end
%     plot(Es, etas);
%     % title("N = 50, Maxwell polynomial basis");
%     str = strcat("(", labels(i), ")");
%     % % attach labels at the tail of each curve
%     % text(Es(end-10*i),etas(end-10*i), str,'fontsize',15);
%     % attach labels at each eta~E curve maximum position
%     text(Es(max_index),etas(max_index), str,'fontsize',15);
%     xlabel('$\varepsilon^* / KT$','Interpreter','latex','fontsize',15);
%     ylabel('$1 - \frac{\lambda_0}{k_{eq}}$','Interpreter','Latex','fontsize',15);
%     legend('$(a): \frac{m}{M}=0.1$', ...
%         '$(b): \frac{m}{M}=0.05$', ...
%         '$(c): \frac{m}{M}=0.02$', ...
%         '$(d): \frac{m}{M}=0.01$', ...
%         'Interpreter','Latex','fontsize',10);
%     hold on;
% end
% savefig("etahat")





%==========================================================================
% % convergence of eta in paper II to eta in paper III when n1/n2 -> 0
% n_ratios = [1, 1/10, 1/100, 1/1000, 1/10000];
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

