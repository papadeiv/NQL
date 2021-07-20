function weights = lagrangian_weights(nodes, I)
    % import and set arbitrary arithmetic precision (vpa) of calculations
    global precision;
    digits(precision);
    % exctract number of nodes and endpoints of the interval of integration
    N = length(nodes);
    a = I(1);
    b = I(2);
    % extract lagrangian basis functions through coefficients matrix
    coefficients = vpa(lagrange(nodes));
    % compute vector of N definite integral in [a,b] of the various monomials of the lagrangian basis functions
    monomial_integral = zeros(N,1);
    monomial_integral = vpa(monomial_integral);
    for k=0:N-1
        % integral of x^k
        monomial_integral(k+1) = vpa((1/(k+1))*(b^(k+1) - a^(k+1)));
    end
    % compute the weigths of the Lagrangian basis
    weights = zeros(N,1);
    for j=1:N
        weights(j) = doubledot(coefficients(j,:),monomial_integral);
    end
end