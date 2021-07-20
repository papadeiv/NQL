function integral = quadrature(kernel, interval, nodes,affine_map_coeff)
    % import and set arbitrary arithmetic precision (vpa) of calculations
    global precision;
    digits(precision);
    alpha = affine_map_coeff(1);
    beta = affine_map_coeff(2);
    N = length(nodes);
    disp('The number of nodes is ')
    disp(N)
    disp('The interpolating Lagrange polynomial and Lagrangian basis function are of degree')
    disp(N-1)
    % compute the weights associated to the nodes
    weights = vpa(lagrangian_weights(nodes, interval));
    % compute the (exact) values of the kernel function in corrispondence of the nodes
    f_values = kernel(alpha*nodes + beta);
    % compute the quadrature integral approximation
    integral = doubledot(weights,f_values);
end