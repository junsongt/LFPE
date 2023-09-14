estars = 1:1:10;
mass_ratio = 0.01;
etas = [];
lambdas = [];
for i = 1:length(estars)
    estar = estars(i);
    result = Lorentz_FPE_with_RX(50,1+1, estar, mass_ratio);
    eta = result.eta;
    etas = [etas eta];
    lambda_0 = result.eigenvalues(1);
    lambdas = [lambdas lambda_0];
end


% eta = g(lambda_0)
plot(lambdas, etas, '-ok');
xlabel('$\lambda_0$','Interpreter','latex','fontsize',20);
ylabel('$\eta$','Interpreter','Latex','fontsize',20);

% % lambda_0 = h(eps)
% plot(estars, lambdas, '-ok');
% xlabel('$\varepsilon^*$','Interpreter','latex','fontsize',20);
% ylabel('$\lambda_0$','Interpreter','Latex','fontsize',20);

% % eps = h(lambda_0)
% plot(lambdas, estars, '-ok');
% xlabel('$\lambda_0$','Interpreter','latex','fontsize',20);
% ylabel('$\varepsilon^*$','Interpreter','Latex','fontsize',20);

%savefig("eta vs lambda_0.fig");