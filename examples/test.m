addpath('../')

close all
clear all
clc

format long

global print_polynomials;
global print_quadrature;
global precision;
global tabulated;

print_polynomials = 1;
print_quadrature = 1;
precision = 32;
tabulated = 1;

a = -pi/4;
b = pi/2;
rule = ["Newton-Cotes";"Gauss-Legendre"];

kernel = @(x) 2*sinh(x.*log(1+x))+(cos(x)).^2;

nodes = [4,12;

for n=1:length(nodes)
    for r=1:2
        main(kernel, rule(r), a, b, nodes(n));
    end
end

