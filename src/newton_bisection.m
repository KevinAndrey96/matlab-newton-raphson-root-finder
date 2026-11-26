function [root, info] = newton_bisection(f, df, a, b, tolerance, max_iterations)
    % NEWTON_BISECTION Safeguarded Newton method on a bracketing interval.
    %
    %   [root, info] = newton_bisection(f, df, a, b, tolerance, max_iterations)
    %
    %   Requires f(a) and f(b) with opposite signs. Each step tries Newton;
    %   when the Newton point leaves the current bracket, or does not halve
    %   the residual fast enough, it falls back to bisection. The bracket
    %   always contains a root, so the method converges for any continuous f
    %   and keeps Newton's quadratic rate near simple roots.
    %
    %   tolerance       width of the final bracket or |f| target (default 1e-12)
    %   max_iterations  iteration budget (default 100)
    %
    %   info fields: converged, iterations, newton_steps, bisection_steps,
    %   residual and bracket ([low, high] at exit).

    if nargin < 5
        tolerance = 1e-12;
    end
    if nargin < 6
        max_iterations = 100;
    end
    if ~isa(f, 'function_handle') || ~isa(df, 'function_handle')
        error('newton:invalidFunction', 'f and df must be function handles');
    end
    if ~isnumeric([a, b]) || numel([a, b]) ~= 2 || any(~isfinite([a, b])) || a >= b
        error('newton:invalidBracket', 'a and b must be finite scalars with a < b');
    end

    fa = f(a);
    fb = f(b);
    if fa == 0
        root = a;
        info = summary(true, 0, 0, 0, 0, [a, a]);

        return
    end
    if fb == 0
        root = b;
        info = summary(true, 0, 0, 0, 0, [b, b]);

        return
    end
    if sign(fa) == sign(fb)
        error('newton:invalidBracket', 'f(a) and f(b) must have opposite signs');
    end

    [low, high] = orient(a, b, fa);
    x = (low + high) / 2;
    fx = f(x);
    previous_step = high - low;
    step = previous_step;
    newton_steps = 0;
    bisection_steps = 0;
    converged = false;
    iteration = 0;

    while iteration < max_iterations
        iteration = iteration + 1;
        slope = df(x);
        newton_x = x - fx / slope;
        leaves_bracket = ~isfinite(newton_x) || (newton_x - low) * (newton_x - high) >= 0;
        too_slow = abs(2 * fx) > abs(previous_step * slope);

        previous_step = step;
        if leaves_bracket || too_slow
            step = (high - low) / 2;
            x = low + step;
            bisection_steps = bisection_steps + 1;
        else
            step = fx / slope;
            x = newton_x;
            newton_steps = newton_steps + 1;
        end

        fx = f(x);
        if fx < 0
            low = x;
        else
            high = x;
        end
        if abs(step) <= tolerance || abs(fx) <= tolerance
            converged = true;
            break
        end
    end

    root = x;
    info = summary(converged, iteration, newton_steps, bisection_steps, abs(fx), ...
                   sort([low, high]));
end

function [low, high] = orient(a, b, fa)
    % Keeps f(low) < 0 < f(high) so the bracket update only checks one sign.
    if fa < 0
        low = a;
        high = b;
    else
        low = b;
        high = a;
    end
end

function info = summary(converged, iterations, newton_steps, bisection_steps, residual, bracket)
    info = struct('converged', converged, 'iterations', iterations, ...
                  'newton_steps', newton_steps, 'bisection_steps', bisection_steps, ...
                  'residual', residual, 'bracket', bracket);
end
