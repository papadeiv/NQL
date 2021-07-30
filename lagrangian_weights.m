function weights = lagrangian_weights(nodes, I)
    global precision;
    % exctract number of nodes and endpoints of the interval of integration
    N = length(nodes);
    a = I(1);
    b = I(2);
    % extract lagrangian basis functions through coefficients matrix
    coefficients = hpf(lagrange(nodes), precision);
    % compute vector of N definite integral in [a,b] of the various monomials of the lagrangian basis functions
    monomial_integral = hpf(zeros(N,1), precision);
    for k=0:N-1
        % integral of x^k
        monomial_integral(k+1) = (hpf('1', precision)/(k+1))*(hpf(sprintf('%.34f', b), precision)^(k+1) - hpf(sprintf('%.34f', a), precision)^(k+1));
    end
    % compute the weigths of the Lagrangian basis
    weights = hpf(zeros(N,1));
    for j=1:N
        weights(j) = hpf(sprintf('%.34f', doubledot(coefficients(j,:),monomial_integral)), precision);
    end
end