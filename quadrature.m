function integral = quadrature(kernel,nodes,affine_map_coeff)
    alpha = affine_map_coeff(1);
    beta = affine_map_coeff(2);
    N = length(nodes);
    disp('The number of nodes is ')
    disp(N)
    disp('The interpolating Lagrange polynomial and Lagrangian basis function are of degree')
    disp(N-1)
    % compute the weights associated to the nodes
    weights = lagrangian_weights(nodes);
    % compute the (exact) values of the kernel function in corrispondence of the nodes
    f_values = kernel(alpha*nodes + beta);
    % compute the quadrature integral approximation
    integral = doubledot(weights,f_values);
end