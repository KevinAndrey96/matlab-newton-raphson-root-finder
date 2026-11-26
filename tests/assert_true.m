function assert_true(condition, label)
    % ASSERT_TRUE Fails with label when condition is false.

    if ~(isscalar(condition) && condition)
        error('tests:false', '%s', label);
    end
end
