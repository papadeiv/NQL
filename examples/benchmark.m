addpath('../')

close all
clear all
clc

format long

global print_polynomials;
global print_quadrature;
global precision;
global tabulated;

print_polynomials = 0;
print_quadrature = 0;
precision = 100;
tabulated = 1;

a = 0;
b = 1;
rule = ["Newton-Cotes";"Gauss-Legendre"];

nodes = [16, 24];
samples = 500;
IEEEaccuracy = eps*ones(samples,1);

for j=1:length(nodes)
    n = nodes(j);
    lambda = linspace(0,2*n,samples);
    apriori = zeros(samples,1);
    aposteriori_exact = zeros(samples,1);
    aposteriori_asymp = zeros(samples,1);
    for k=1:samples
        tic
        % define the integrand function (Muntz monomial of non-integer degree lambda(r))
        kernel = @(x) x.^lambda(k);
        % compute a-priori estimate of Muntz monomial
        I_h = main(kernel, rule(2), a, b, n); % approximate integral using the specified quadrature rule
        I = 1/(lambda(k)+1); % exact definite intergal of Muntz monomial in [0,1]
        apriori(k) = abs(I - I_h);
        % compute a-posteriori exact and asymptotic estimates
        z = 2*lambda(k);
        w = 2*n - lambda(k);
        t = 2 + 2*lambda(k);
        den = 2*n + lambda(k);
        common_factor = -2^(1+lambda(k))*2*sin(pi*lambda(k));
        aposteriori_exact(k) = common_factor*2^(-lambda(k))*lambda(k)*(beta(z,w)/(den) - beta(z,2 + w)/(2 + den));
        aposteriori_asymp(k) = common_factor*2^(-lambda(k))*(1+2*n)^(-t)*gamma(t);
        toc
    end
    text = ['n = ', num2str(n)];
    figure
    semilogy(lambda,IEEEaccuracy, 'color', 'k','LineWidth', 2.5, 'DisplayName', 'IEEE double precision');
    hold on
    semilogy(lambda,apriori, 'color', 'r','LineWidth', 2.5, 'DisplayName', 'A-priori estimate');
    semilogy(lambda,aposteriori_exact, 'color', 'b', 'LineWidth', 2.5, 'DisplayName', 'A-posteriori exact estimate');
    semilogy(lambda,aposteriori_asymp, '--', 'color', '#7E2F8E', 'LineWidth', 2.5, 'DisplayName', 'A-posteriori asymptotic estimate');
    leg = legend('Interpreter', 'latex', 'Interpreter', 'latex', 'Interpreter', 'latex');
    set(leg, 'Location', 'best', 'FontSize', 18);
    xlabel('\lambda')
    title(text)
    grid on
    hold off
end

