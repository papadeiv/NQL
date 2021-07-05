function nodes = legendre(N)
    global print_flag;
    % Import colors from custom-defined colormap
    color = cmap(N+1,1);
    % Define coefficient matrix
    coefficients = zeros(N+1,N+1);
    % Assign P0(x)
    coefficients(1,1) = 1;
    % Assign P1(x)
    coefficients(2,1) = 0;
    coefficients(2,2) = 1;
    % Construct the coefficient matrix
    for n=3:N+1
        A = (2*n+1)/(n+1);
        B = 0;
        C = n/(n+1);
        coefficients(n,2:n) = A*coefficients(n-1,1:n-1);
        coefficients(n,:) = coefficients(n,:) - C*coefficients(n-2,:);
    end
    % Compute the values of monomials over [-1,1]
    x = linspace(-1,1,100000);
    monomials = zeros(length(x),N+1);
    for n=0:N
        monomials(:,n+1) = x.^n;
    end
    if print_flag==1
        figure('units','normalized','outerposition',[0 0 1 1])
    end
    for n=0:N
        % Compute the values of Legendre's n-th polynomial over [-1,1]
        polynomials = zeros(length(x),1);
        for j=1:length(x)
            polynomials(j) = doubledot(coefficients(n+1,:),monomials(j,:));
        end
        if print_flag==1
            % Plot Legendre's n-th polynomial
            hold on
            grid on
            p(n+1) = plot(x,polynomials,'LineWidth',3, 'color', color(n+1,:), 'DisplayName',sprintf('$ P_{%d}(x) $', n));
        end
        % Find and extract indeces of zero-crossing values of the polynomial
        % array i.e. the contiguous cells of the array that change sign (as per Rolle's theorem)
        counter = 0;
        for k=2:length(polynomials)
            if polynomials(k-1)*polynomials(k) <= 0
                counter = counter + 1;
                indeces(counter,:) = [k-1,k];
            end
        end
        nodes = zeros(counter,1);
        % Find the values in x associated to the zero-crossing cells
        % and interpolate linearly to compute the associated nodes
        for j=1:counter
            nodes(j) = x(indeces(j,1)) + (polynomials(indeces(j,1))*(x(indeces(j,1))-x(indeces(j,2))))/(polynomials(indeces(j,1))-polynomials(indeces(j,2)));         
        end
        if print_flag==1
            % Plot the nodes of the Legendre's n-th polynomial
            plot(nodes,zeros(length(nodes)),'o','MarkerSize',10, 'MarkerEdgeColor','k','MarkerFaceColor', color(n+1,:));
        end
    end
    if print_flag==1
        leg = legend(p,'Interpreter', 'latex','Orientation','horizontal','NumColumns',4);
        set(leg, 'Location', 'bestoutside', 'FontSize', 15)
    end
end