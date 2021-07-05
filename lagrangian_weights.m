function weights = lagrangian_weights(nodes)
    % exctract number of nodes and endpoints of the interval of integration
    N = length(nodes);
    a = nodes(1);
    b = nodes(N);
    % extract lagrangian basis functions through coefficients matrix
    coefficients = lagrange(nodes);
    % compute vector of N definite integral in [a,b] of the various monomials of the lagrangian basis functions
    monomial_integral = zeros(N,1);
    for k=0:N-1
        % integral of x^k
        monomial_integral(k+1) = (1/(k+1))*(b^(k+1) - a^(k+1));
    end
    % compute the weigths of the Newton-Cotes quadrature rule
    weights = zeros(N,1);
    for j=1:N
        weights(j) = doubledot(coefficients(j,:),monomial_integral);
    end 
end