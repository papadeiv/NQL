function basis = lagrange(nodes)
    global precision;
    % extract number of nodes
    N = length(nodes);
    % Import colors from custom-defined colormap
    color = cmap(N,1);
    % construct all the lagrangian monomial iteratively and extract the lagrangian basis of the (N-1)-th order
    basis = hpf(zeros(N,N));
    for j=1:N
        % of all the nodes in [a,b] remove the one associated to the j-th basis function (in which its value has to be 1)
        lagrangian_nodes = nodes([1:j-1,j+1:N]);
        % on the j-th row store the N coefficients of the j-th lagrangian polynomial of degree (N-1)
        basis(j,:) = hpf(extract_basis(lagrangian_nodes,nodes(j)), precision); % the first row obviously stores the coefficient of the power x^0
    end
    
    global print_polynomials;
    global print_quadrature;
    
    % check consinstency in the print options selection
    if print_quadrature==1 && print_polynomials==0
        error('In order to print the Lagrangian interpolating polynomial (print_quadrature=1) you must also enable the option to print Lagrangian basis functions (print_polynomials=1).');
    end
    
    if print_polynomials==1
        % initialise the monomial matrix in [a,b]
        x = linspace(nodes(1), nodes(N), 1000);
        monomials = zeros(length(x),N);
        % fill-in the rows of the matrix by the value of the monomial power
        for n=1:N
            % along the n-th column store the value of the vector x = (0.000, 0.100, 0.200, ...) to the power (n-1) to compute x^0, x^1, x^2 ...
            monomials(:,n) = x.^(n-1);
        end
        
        if print_quadrature==1
            global polynomials;
        end
        
        polynomials = zeros(length(x),N);
        
        figure('units','normalized','outerposition',[0 0 1 1])
        % to represent the polynomial basis compute the scalar product between the coefficients in the basis matrix and the values of (x^0, x^1, x^2, ...) associated to each x in [a,b]
        for n=1:N
            for j=1:length(x)
                polynomials(j,n) = doubledot(basis(n,:), monomials(j,:));
            end
            hold on
            grid on
            p(n)=plot(x, polynomials(:,n), 'LineWidth', 3, 'color', color(n,:), 'DisplayName',sprintf('$ l_{%d}(x) $', n-1));
        end
        for n=1:N
            plot(nodes(n), 0, 'o','MarkerSize',10, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', color(n,:));
        end

        leg = legend(p,'Interpreter', 'latex','Orientation','horizontal','NumColumns',4);
        set(leg, 'Location', 'bestoutside', 'FontSize', 15);

        str = sprintf('Lagrangian basis functions: $l_{j}(x)$ are $P_{%d}$ polynomials',N-1);
        title(str, 'interpreter', 'latex', 'FontSize', 14)
    end
end