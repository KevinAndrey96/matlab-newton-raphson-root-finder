% NEWTON_DEMO Shows Newton's method on the original 2019 exercises and on a
% repeated root, printing every iterate and the observed convergence order.
%
% Run from the repository root:
%   octave-cli examples/newton_demo.m
%   matlab -batch "run('examples/newton_demo.m')"

here = fileparts(mfilename('fullpath'));
addpath(here, fullfile(here, '..', 'src'));

print_newton_run('sin(x) from x0 = 3 (original exercise)', @sin, 3, struct(), pi);

polynomial = @(x) x.^2.01 - 13 * x + 30;
[low_root, ~] = newton_raphson(polynomial, 2);
print_newton_run('x^2.01 - 13x + 30 from x0 = 2', polynomial, 2, struct(), low_root);

double_root = @(x) (x - 1).^2;
plain = struct('derivative', @(x) 2 * (x - 1), 'tolerance', 1e-14);
print_newton_run('(x - 1)^2 without multiplicity', double_root, 3, plain, 1);
corrected = plain;
corrected.multiplicity = 2;
print_newton_run('(x - 1)^2 with multiplicity 2', double_root, 3, corrected, 1);

[root, info] = newton_bisection(@atan, @(x) 1 ./ (1 + x.^2), -10, 3);
fprintf('\natan(x) on [-10, 3] with the safeguarded method: root %.3e\n', root);
fprintf('  %d Newton steps, %d bisection steps (plain Newton from 3 diverges)\n', ...
        info.newton_steps, info.bisection_steps);
