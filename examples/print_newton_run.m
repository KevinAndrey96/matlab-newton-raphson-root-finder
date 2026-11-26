function print_newton_run(title, f, x0, options, exact)
    % PRINT_NEWTON_RUN Prints every Newton iterate with its error and order.

    [root, info] = newton_raphson(f, x0, options);
    orders = [NaN, NaN, convergence_order(info.history, exact)];
    fprintf('\n%s\n', title);
    fprintf('%4s  %22s  %12s  %7s\n', 'k', 'x(k)', '|x(k) - r|', 'order');
    for k = 1:numel(info.history)
        order = NaN;
        if k <= numel(orders)
            order = orders(k);
        end
        fprintf('%4d  %22.15f  %12.3e  %7.3f\n', k - 1, info.history(k), ...
                abs(info.history(k) - exact), order);
    end
    fprintf('root %.15f after %d iterations, |f(root)| = %.3e\n', root, info.iterations, ...
            info.residual);
end
