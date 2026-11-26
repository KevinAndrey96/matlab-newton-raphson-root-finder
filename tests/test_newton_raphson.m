function tests = test_newton_raphson()
    tests = {@finds_pi_from_the_original_start, @uses_an_analytic_derivative, ...
             @converges_quadratically_on_simple_roots, @is_linear_on_double_roots, ...
             @multiplicity_restores_quadratic_rate, @solves_the_original_polynomial, ...
             @returns_immediately_on_a_root, @reports_when_the_budget_runs_out, ...
             @rejects_zero_derivatives, @rejects_diverging_iterations, @validates_input};
end

function finds_pi_from_the_original_start()
    [root, info] = newton_raphson(@sin, 3);
    assert_close(root, pi, 1e-12, 'root of sin near 3');
    assert_true(info.converged && info.iterations <= 5, 'converges in a few steps');
    assert_true(info.residual <= 1e-10, 'residual below tolerance');
end

function uses_an_analytic_derivative()
    options = struct('derivative', @(x) 2 * x, 'tolerance', 1e-15);
    [root, info] = newton_raphson(@(x) x.^2 - 2, 1, options);
    assert_close(root, sqrt(2), 1e-12, 'sqrt(2)');
    assert_close(info.history(2), 1.5, 0, 'first Newton step');
end

function converges_quadratically_on_simple_roots()
    options = struct('derivative', @(x) 3 * x.^2 - 2, 'tolerance', 1e-15);
    [root, info] = newton_raphson(@(x) x.^3 - 2 * x - 5, 2, options);
    orders = convergence_order(info.history, root);
    assert_true(~isempty(orders) && abs(orders(end) - 2) < 0.2, ...
                sprintf('last order %s should be close to 2', mat2str(orders, 3)));
end

function is_linear_on_double_roots()
    options = struct('derivative', @(x) 2 * (x - 1), 'tolerance', 1e-14);
    [root, info] = newton_raphson(@(x) (x - 1).^2, 3, options);
    assert_close(root, 1, 1e-6, 'double root');
    assert_true(info.iterations > 15, 'needs many steps without multiplicity');
    orders = convergence_order(info.history, 1);
    assert_true(abs(orders(end) - 1) < 0.05, 'linear convergence');
end

function multiplicity_restores_quadratic_rate()
    options = struct('derivative', @(x) 3 * (x - 2).^2, 'multiplicity', 3);
    [root, info] = newton_raphson(@(x) (x - 2).^3, 5, options);
    assert_close(root, 2, 1e-12, 'triple root');
    assert_true(info.iterations <= 2, 'exact after one corrected step');
end

function solves_the_original_polynomial()
    f = @(x) x.^2.01 - 13 * x + 30;
    [low, low_info] = newton_raphson(f, 2);
    [high, high_info] = newton_raphson(f, 11);
    assert_true(low_info.converged && high_info.converged, 'both roots converge');
    assert_true(abs(f(low)) <= 1e-10 && abs(f(high)) <= 1e-10, 'residuals');
    assert_true(low < 4 && high > 9, 'two distinct roots');
end

function returns_immediately_on_a_root()
    [root, info] = newton_raphson(@sin, 0);
    assert_close(root, 0, 0, 'start at a root');
    assert_true(info.converged && info.iterations == 0, 'no steps');
end

function reports_when_the_budget_runs_out()
    options = struct('max_iterations', 3, 'tolerance', 1e-30, 'step_tolerance', 1e-30);
    [~, info] = newton_raphson(@(x) x.^2 - 2, 100, options);
    assert_true(~info.converged && info.iterations == 3, 'stops at the budget');
    assert_true(numel(info.history) == 4, 'history keeps x0 and every step');
end

function rejects_zero_derivatives()
    options = struct('derivative', @(x) 0 * x);
    assert_error_id(@() newton_raphson(@(x) x.^2 + 1, 0.5, options), 'newton:zeroDerivative');
end

function rejects_diverging_iterations()
    cube_root = @(x) sign(x) .* abs(x).^(1 / 3);
    options = struct('derivative', @(x) abs(x).^(-2 / 3) / 3, 'max_iterations', 5000);
    assert_error_id(@() newton_raphson(cube_root, 1, options), 'newton:nonFinite');
end

function validates_input()
    assert_error_id(@() newton_raphson('sin', 3), 'newton:invalidFunction');
    assert_error_id(@() newton_raphson(@sin, Inf), 'newton:invalidStart');
    assert_error_id(@() newton_raphson(@sin, 3, struct('speed', 2)), 'newton:invalidOptions');
    assert_error_id(@() newton_raphson(@sin, 3, struct('tolerance', -1)), 'newton:invalidOptions');
    assert_error_id(@() newton_raphson(@sin, 3, struct('max_iterations', 2.5)), ...
                    'newton:invalidOptions');
end
