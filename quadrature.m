function integral = quadrature(kernel, interval, nodes, affine_map_coeff)
    global precision;
    alpha = affine_map_coeff(1);
    beta = affine_map_coeff(2);
    N = length(nodes);
    % compute the weights associated to the nodes
    weights = hpf(lagrangian_weights(nodes, interval), precision);
    % compute the (exact) values of the kernel function in corrispondence of the nodes
    f_values = kernel(alpha*nodes + beta);
    % compute the quadrature integral approximation
    integral = hpf(sprintf('%.34f', doubledot(weights,f_values)), precision);
    % print the Lagrangian interpolating polynomial if option is selected
    global print_quadrature;
    if print_quadrature==1
        x = linspace(nodes(1), nodes(N), 1000);
        global polynomials;
        for n=1:N
            polynomials(:,n) = polynomials(:,n).*f_values(n);
        end
        lagrangian_interpolating_polynomial = sum(transpose(polynomials));
        
        figure('units','normalized','outerposition',[0 0 1 1])
        hold on
        grid on
        func = plot(x, kernel(alpha*x + beta), '--', 'color', 'r', 'LineWidth', 2);
        interp = plot(x,lagrangian_interpolating_polynomial, 'color', 'k', 'LineWidth', 2.5);
        plot(nodes, f_values, 'o', 'color', 'k', 'MarkerFaceColor', 'k', 'MarkerSize', 8);
        legend([interp, func], {'$ L_n (x) $', '$ f(x) $'}, 'Interpreter', 'latex', 'Location', 'best', 'FontSize', 15);
    end
end