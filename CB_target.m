function err = CB_target(K, L, Y, x, tau)
% целевая функция - минимум суммы квадратов со скользящим окном
K = K(tau:end); L = L(tau:end); Y = Y(tau:end);
Y_fit = x(1) .* K .^ x(2) .* L .^ x(3);
err = sum((Y_fit - Y) .^ 2);