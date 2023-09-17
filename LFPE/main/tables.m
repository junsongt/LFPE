% Author: Junsong Tang
% This script could generate:
% (1) Table with convergence of eigenvalues of LFP opearator: lambda_i 
% of selected order i with increasing N;

% (2) Table of convergence of eta with different mass ratio setting
% in increasing number of Sonine polynomials;

% (3) Table of lambda0 with different mass ratios and different E/KT;

% (4) Table of lambda0 with decreasing mass ratios with different fixed E/KT

% plot the lower order eigenfunctions: psi_0, psi_1, psi_2,...psi_5

% Note that if users want to choose one of the features above, make sure to
% comment out the other blocks


close; clear; clc;


% Ns = [5,10,15,20,30,40,50];
% Ns = [1,2,4,6,8,10];
% Ns = [1,2,3,4,5,6,8,10];

%==========================================================================
% (1) Table of convergence of eigenvalues of LFP operator(w/wo R(x))
estar = 2;
mass_ratio = 1/10;

Ns = [5,10,15,20,30,40,50];
% selected order j of lambda_j: [lam_0, lam_1, lam_2, ...lam_20]
selected_index = [0,1,2,3,4,5,10,15,20];

% get eigenvalues for different N to check convergence of eigenvalues
table = [];
for i = 1:length(Ns)
    % result = LFPE(Ns(i)+1,1+1);
    result = LFPE_Rx(Ns(i)+1,1+1, estar, mass_ratio);
    lambdas = result.eigenvalues_Rx_ps;
    selected_lambdas = [];
    for j = 1:length(selected_index)
        lambda_idx = selected_index(j)+1;
        if lambda_idx > length(lambdas)
            lambda = NaN;
        else
            lambda = lambdas(lambda_idx);
        end
        selected_lambdas = [selected_lambdas lambda];
    end
    table(i,:) = [Ns(i) selected_lambdas];
end
% input converted to be latex form
input.data = table;
% input.dataFormat = {'%.4f'};
input.dataFormat = {'%.4g'};



%==========================================================================
% (2) Table of convergence of eta with different mass ratio setting
% in increasing number of Sonine polynomials
% mass_ratios = [1, 0.1, 0.01, 0.001];
mass_ratios = [2, 1, 0.8, 0.5, 0.1, 0.05];
% mass_ratios = [0.6, 0.3, 0.15, 0.05, 0.02, 0.01];
table = [];
E = 15;
rn = 1;
% delta := (dR/dE)^2
delta = 1;
for i = 1:length(Ns)
    N = Ns(i);
    etas = [];
    for j = 1:length(mass_ratios)
        rm = mass_ratios(j);
        result = NEQ3(rm, E, rn, delta, N);
        eta = result.eta;
        etas = [etas eta];
    end
    etas = 1000 * etas;
    table(i,:) = [N etas];
end
% input converted to be latex form
input.data = table;
% input.dataFormat = {'%.4f'};
input.dataFormat = {'%.4g'};


%==========================================================================
% (3) Table of lambda0 with different mass ratios and different E/KT
N = 50;
mass_ratios = [1, 0.5, 0.1, 0.05, 0.01, 0.005, 0.001];
table = [];
Es = [1,2,3,4,6,8,10];
for i = 1:length(Es)
    lams = [];
    E = Es(i);
    for j = 1:length(mass_ratios)
        rm = mass_ratios(j);
        result = LFPE_Rx(N,1,E,rm);
        lam0 = result.lambda0_R;
        lams(j) = lam0;
    end
    table(i, :) = [E, lams];
end
% input converted to be latex form
input.data = table;
% input.dataFormat = {'%.4f'};
input.dataFormat = {'%.4g'};



%==========================================================================
% (4) Table of lambda0 with decreasing mass ratios with different fixed E/KT
% this table we could see the sign change in lambda_0 for some critical
% mass ratio
N = 50;
table = [];
Es = [2,3,4,6,8,10];
mass_ratios = [1, 0.5, 0.2, 0.15, 0.1, 0.05, 0.02, 0.01];
for i = 1:length(mass_ratios)
    lams = [];
    rm = mass_ratios(i);
    for j = 1:length(Es)
        E = Es(j);
        result = LFPE_Rx(N,1,E,rm);
        lam0 = result.lambda0_R;
        lams(j) = lam0;
    end
    table(i, :) = [rm, lams];
end

% table = array2table(table, "VariableNames", ["$N$","$\lambda_1$", ...
%     "$\lambda_2$", "$\lambda_3$", "$\lambda_4$", "$\lambda_5$", ...
%     "$\lambda_10$", "$\lambda_15$", "$\lambda_20$"]);

% input converted to be latex form
input.data = table;
% input.dataFormat = {'%.4f'};
input.dataFormat = {'%.4g'};



%==========================================================================
% (5) generate lower-order eigenfunctions given number of basis = 50
n = 50;
for order = 0:5
    %result = Lorentz_FPE(n+1,order+1);
    % result = Lorentz_FPE_with_RX(n+1, order+1, estar, mass_ratio);
    result = LFPE_Rx(n+1, order+1, estar, mass_ratio);
    eigfcn = result.eigenfunction_Rx;
    pt = result.quadpt;
    plot_eigenfunction(pt, eigfcn, order+1);
end


% Helper function to plot the npoly-th eigenfunctions
function plot_eigenfunction(pt,eigfcn, npoly)
fac=eigfcn(npoly)/abs(eigfcn(npoly));
plot(pt,eigfcn,'-ok','linewidth',1.2,'markersize',5,'markerfacecolor','k')
axis([0 5 -.5 .5])
set(gca,'FontSize',24)
set(gca,'Ytick',[-.4:.2:.4],'linewidth',1.6)
set(gca,'Xtick',[0:1:5],'linewidth',1.6)
str2=num2str(npoly-1);
str1={'$n  =  $'};
str=strcat(str1, str2);
text(3,-.4,0,str,'Interpreter','latex','fontsize',26)
xlabel('$x$','Interpreter','latex','fontsize',30)
ylabel('$\psi_n(x)$','Interpreter','Latex','fontsize',30)
savefig("psi"+str2+".fig");
end
