function  integral = main(kernel, rule, a, b, N)
    global precision;
    disp('The number of nodes is ')
    disp(N)
    disp('The Lagrange interpolating polynomial and Lagrangian basis functions have degree')
    disp(N-1)
    % define vector containing coefficients of affine map [-1,1] -> [a,b] and initialise for no map case (i.e. Newton-Cotes formulae)
    coefficients = [1, 0];
    if rule=="Newton-Cotes"
        % define the jacobian associated to the affine map and initialise for no map case
        jacobian = 1;
        % build equispaced nodes in [a,b]
        nodes = linspace(hpf(a, precision), hpf(b, precision), N);
        % define the integration interval
        I = [a, b];
        disp('Newton-Cotes quadrature formula on equispaced nodes will be used which has an order of precision')
        disp(N-1)
    else
        % compute the coefficients of the linear map associated to the interval [a,b] for the G-L rule
        coefficients(1) = (hpf('1', precision)/2)*(b-a);
        coefficients(2) = (hpf('1', precision)/2)*(a+b);
        % compute the jacobian of the affine map associated to the interval [a,b] for the G-L rule
        jacobian = (hpf(sprintf('%.34f', b), precision) - hpf(sprintf('%.34f', a), precision))/2;
        % numerically compute G-L nodes approximation
        nodes = legendre(N);
        % define the integration interval
        I = [-1, 1];
        disp('Gauss-Legendre quadrature formula will be used which has the roots of Legendre polynomials as nodes and whose order of precision is')
        disp(2*N-1)
    end
    integral = jacobian*hpf(sprintf('%.34f', quadrature(kernel, I, nodes, coefficients)), precision);
    sprintf('The approximated integral is I_h = %f', integral)
end