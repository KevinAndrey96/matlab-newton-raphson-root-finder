function orders = convergence_order(history, root)
    % CONVERGENCE_ORDER Estimates the order p of an iteration from its history.
    %
    %   orders = convergence_order(history, root)
    %
    %   With errors e(k) = |x(k) - root|, each estimate is
    %   p(k) = log(e(k+1) / e(k)) / log(e(k) / e(k-1)).
    %   Newton gives values close to 2 near a simple root and close to 1 near
    %   a repeated root. Steps whose error already reached zero are skipped.

    errors = abs(history(:)' - root);
    errors = errors(errors > 0);
    if numel(errors) < 3
        orders = [];

        return
    end
    orders = log(errors(3:end) ./ errors(2:end - 1)) ./ log(errors(2:end - 1) ./ errors(1:end - 2));
end
