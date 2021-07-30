function polynomial = extract_basis(nodes,pivot)
global precision;
    % construct the lagrangian polynomial of degree N-1 iteratively by ascending degree order
    N = length(nodes)+1;
    if N < 2
        error('Lagrange basis functions have to be associated to at least 2 nodes!')
    end
    % initialise the coefficients matrix
    coefficients = hpf(zeros(N-1,N));
    % along the n-th row store the values of the N coefficients of the n-th degree polynomial
    coefficients(1,1) = hpf(sprintf('%.34f', -nodes(1)), precision);
    coefficients(1,2) = hpf(sprintf('%.34f', 1), precision);
    % compute the numerator of the (N-1)-th degree lagrangian polynomial iteratevely by ascending degree order
    for n=2:(N-1)
        coefficients(n,2:n+1) = coefficients(n-1,1:n);
        coefficients(n,:) = coefficients(n,:) - hpf(sprintf('%.34f', nodes(n)), precision).*coefficients(n-1,:);
    end
    % compute the denominator of the (N-1)-degree lagrangian polynomial by multiplying all the denomiantor associated to each node
    denominator = hpf('1', precision);
    for n=1:(N-1)
        denominator = denominator*(hpf(sprintf('%.34f', pivot), precision) - hpf(sprintf('%.34f', nodes(n)), precision));
    end
    % return the last row of the (scaled) matrix as Lagrange's (N-1)-th degree polynomial associated to the input pivot and nodes selection
    polynomial = (hpf('1', precision)/denominator)*coefficients(N-1,:);
end