function err = CB_target(K, L, Y, x, w)
% целевая функция - минимум суммы квадратов со скользящим окном
K = K(1:w); L = L(1:w); Y = Y(1:w);
Y_mod = x(1) .* K .^ x(2) .* L .^ x(3);
err = sum((Y_mod - Y) .^ 2);