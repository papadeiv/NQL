addpath('../')

close all
clear all
clc

format longG

global print_polynomials;
global print_quadrature;
global precision;
global tabulated;

print_polynomials = 0;
print_quadrature = 0;
precision = 16;
tabulated = 1;

a = 0;
b = 1;
rule = ["Newton-Cotes";"Gauss-Legendre"];

nodes = [24];
samples = 1000;

IEEEaccuracy = eps*ones(samples,1);

for j=1:length(nodes)
    n = nodes(j);
    lambda = linspace(hpf('0', precision),hpf('40', precision),samples);
    exact_error = hpf(zeros(samples,1), precision);
    aposteriori_exact = hpf(zeros(samples,1), precision);
    aposteriori_asymp = hpf(zeros(samples,1), precision);
    for k=1:samples
        tic
        % define the integrand function (Muntz monomial of non-integer degree lambda(r))
        kernel = @(x) x.^lambda(k);
        % compute a-priori estimate of Muntz monomial
        I_h = main(kernel, rule(2), a, b, n); % approximate integral using the specified quadrature rule
        I = 1/(lambda(k)+1); % exact definite intergal of Muntz monomial in [0,1]
        exact_error(k) = abs(I- I_h);
        % compute a-posteriori exact and asymptotic estimates
        z = 2*double(lambda(k));
        w = 2*n - double(lambda(k));
        t = 2 + 2*double(lambda(k));
        den = 2*n + lambda(k);
        beta1 = (hpf(gamma(z),precision)*hpf(gamma(w),precision))/hpf(gamma(z+w), precision);
        beta2 = (hpf(gamma(z),precision)*hpf(gamma(2+w),precision))/hpf(gamma(z+2+w), precision);
        common_factor = -2^(1+lambda(k))*2*hpf(sin(hpf('pi', precision)*lambda(k)), precision);
        aposteriori_exact(k) = abs(common_factor*hpf(2^(double(-lambda(k))), precision)*lambda(k)*(hpf(beta1/(den), precision) - hpf(beta2/(2 + den), precision)));
        aposteriori_asymp(k) = abs(common_factor*hpf(2^(double(-lambda(k))), precision)*hpf((1+2*n)^(-t), precision)*hpf(gamma(t), precision));
        toc
    end
    text = ['n = ', num2str(n)];
    figure
    semilogy(lambda,IEEEaccuracy, 'color', 'k','LineWidth', 2.5, 'DisplayName', 'IEEE double precision');
    hold on
    semilogy(lambda,exact_error, 'color', 'r','LineWidth', 2.5, 'DisplayName', 'Exact error');
    semilogy(lambda,aposteriori_exact, 'color', 'b', 'LineWidth', 2.5, 'DisplayName', 'A-posteriori exact estimate');
    semilogy(lambda,aposteriori_asymp, '--', 'color', '#7E2F8E', 'LineWidth', 2.5, 'DisplayName', 'A-posteriori asymptotic estimate');
    leg = legend('Interpreter', 'latex', 'Interpreter', 'latex', 'Interpreter', 'latex');
    set(leg, 'Location', 'best', 'FontSize', 18);
    xlabel('\lambda')
    title(text)
    grid on
    hold off
end


