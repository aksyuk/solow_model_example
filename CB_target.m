function err = CB_target(K, L, Y, x, w)
% целевая функция - минимум суммы квадратов со скользящим окном
K = K(w:end); L = L(w:end); Y = Y(w:end);
Y_fit = x(1) .* K .^ x(2) .* L .^ x(3);
err = sum((Y_fit - Y) .^ 2);