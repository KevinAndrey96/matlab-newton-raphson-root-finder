function [root, info] = newton_raphson(f, x0, options)
    % NEWTON_RAPHSON Finds a root of f starting from x0 with Newton's method.
    %
    %   [root, info] = newton_raphson(f, x0)
    %   [root, info] = newton_raphson(f, x0, options)
    %
    %   Iterates x(k+1) = x(k) - m * f(x(k)) / f'(x(k)) until |f(x)| or the
    %   step size falls below its tolerance.
    %
    %   options fields (all optional):
    %     derivative      handle to f'; a central difference is used when empty
    %     tolerance       residual tolerance on |f(x)| (default 1e-10)
    %     step_tolerance  relative tolerance on |x(k+1) - x(k)| (default 1e-14)
    %     max_iterations  iteration budget (default 50)
    %     multiplicity    known root multiplicity m; m > 1 restores quadratic
    %                     convergence on repeated roots (default 1)
    %
    %   info fields:
    %     converged   true when a tolerance was met within the budget
    %     iterations  Newton steps taken
    %     residual    |f(root)|
    %     history     every iterate, starting with x0
    %
    %   Errors: newton:invalidFunction, newton:invalidStart,
    %   newton:invalidOptions, newton:zeroDerivative, newton:nonFinite

    if nargin < 3
        options = struct();
    end
    if ~isa(f, 'function_handle')
        error('newton:invalidFunction', 'f must be a function handle');
    end
    if ~isnumeric(x0) || ~isscalar(x0) || ~isreal(x0) || ~isfinite(x0)
        error('newton:invalidStart', 'x0 must be a finite real scalar');
    end
    settings = resolve_options(options);

    x = x0;
    fx = f(x);
    history = zeros(1, settings.max_iterations + 1);
    history(1) = x;
    iterations = 0;
    converged = abs(fx) <= settings.tolerance;

    while ~converged && iterations < settings.max_iterations
        slope = derivative_at(f, settings.derivative, x);
        if slope == 0
            error('newton:zeroDerivative', 'derivative vanished at x = %.17g', x);
        end
        step = settings.multiplicity * fx / slope;
        x = x - step;
        fx = f(x);
        if ~isfinite(x) || ~isfinite(fx)
            error('newton:nonFinite', 'iteration diverged after %d steps', iterations + 1);
        end
        iterations = iterations + 1;
        history(iterations + 1) = x;
        converged = abs(fx) <= settings.tolerance || ...
                    abs(step) <= settings.step_tolerance * max(1, abs(x));
    end

    root = x;
    info = struct('converged', converged, 'iterations', iterations, ...
                  'residual', abs(fx), 'history', history(1:iterations + 1));
end

function settings = resolve_options(options)
    if ~isstruct(options) || ~isscalar(options)
        error('newton:invalidOptions', 'options must be a scalar struct');
    end
    settings = struct('derivative', [], 'tolerance', 1e-10, 'step_tolerance', 1e-14, ...
                      'max_iterations', 50, 'multiplicity', 1);
    names = fieldnames(options);
    for k = 1:numel(names)
        if ~isfield(settings, names{k})
            error('newton:invalidOptions', 'unknown option %s', names{k});
        end
        settings.(names{k}) = options.(names{k});
    end
    if ~isempty(settings.derivative) && ~isa(settings.derivative, 'function_handle')
        error('newton:invalidOptions', 'derivative must be a function handle');
    end
    positive = [settings.tolerance, settings.step_tolerance, settings.max_iterations, ...
                settings.multiplicity];
    if ~isnumeric(positive) || numel(positive) ~= 4 || any(~(positive > 0)) || ...
       mod(settings.max_iterations, 1) ~= 0
        error('newton:invalidOptions', ...
              'tolerances and multiplicity must be positive and max_iterations a positive integer');
    end
end

function slope = derivative_at(f, derivative, x)
    if ~isempty(derivative)
        slope = derivative(x);

        return
    end
    h = eps^(1 / 3) * max(1, abs(x));
    slope = (f(x + h) - f(x - h)) / (2 * h);
end
