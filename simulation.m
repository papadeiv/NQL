close all
clear all
clc

global print_flag;

a = 0;
b = 1;
N = 10;
rule = ["Newton-Cotes";"Gauss-Legendre"];

integrand = @(x) x.^0;

computed_integral = zeros(N-1,1);

for n=2:N
    print_flag = (N+1) - n;
    computed_integral(n-1) = main(integrand, rule(2), a, b, n);
end
figure
plot(computed_integral, '-o', 'LineWidth', 2.5, 'color', 'r')

