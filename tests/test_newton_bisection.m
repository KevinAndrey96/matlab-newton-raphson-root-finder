function tests = test_newton_bisection()
    tests = {@converges_where_plain_newton_diverges, @matches_newton_on_well_behaved_roots, ...
             @accepts_decreasing_functions, @returns_endpoints_that_are_roots, ...
             @keeps_the_root_inside_the_bracket, @rejects_invalid_brackets};
end

function converges_where_plain_newton_diverges()
    [root, info] = newton_bisection(@atan, @(x) 1 ./ (1 + x.^2), -10, 3);
    assert_close(root, 0, 1e-12, 'root of atan');
    assert_true(info.converged && info.bisection_steps > 0, 'needs bisection safeguards');
end

function matches_newton_on_well_behaved_roots()
    [root, info] = newton_bisection(@(x) x.^3 - 2 * x - 5, @(x) 3 * x.^2 - 2, 2, 3);
    assert_close(root, 2.0945514815423265, 1e-12, 'Wallis cubic');
    assert_true(info.newton_steps > info.bisection_steps, 'mostly Newton steps');
end

function accepts_decreasing_functions()
    [root, info] = newton_bisection(@(x) cos(x) - x, @(x) -sin(x) - 1, 0, 1);
    assert_close(root, 0.7390851332151607, 1e-12, 'Dottie number');
    assert_true(info.converged, 'converged');
end

function returns_endpoints_that_are_roots()
    [root, info] = newton_bisection(@(x) x - 1, @(x) 1, 1, 2);
    assert_close(root, 1, 0, 'left endpoint');
    assert_true(info.iterations == 0, 'no iterations');
    root = newton_bisection(@(x) x - 2, @(x) 1, 1, 2);
    assert_close(root, 2, 0, 'right endpoint');
end

function keeps_the_root_inside_the_bracket()
    [root, info] = newton_bisection(@(x) x.^3, @(x) 3 * x.^2, -1, 2, 1e-10);
    assert_true(abs(root) < 1e-3, 'triple root at zero');
    assert_true(info.bracket(1) <= root && root <= info.bracket(2), 'bracket contains root');
end

function rejects_invalid_brackets()
    assert_error_id(@() newton_bisection(@sin, @cos, 1, 2), 'newton:invalidBracket');
    assert_error_id(@() newton_bisection(@sin, @cos, 2, 1), 'newton:invalidBracket');
    assert_error_id(@() newton_bisection(@sin, @cos, -Inf, 1), 'newton:invalidBracket');
    assert_error_id(@() newton_bisection(@sin, 'cos', -1, 1), 'newton:invalidFunction');
end
