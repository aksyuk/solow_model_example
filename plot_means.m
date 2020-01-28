function [m_all, m_end] = plot_means(x0, y0, str_title, str_legend_1, ...
    file_name)

clf
hold on
n = length(x0); x1 = x0; x1(1:19) = nan; 
m_end = mean(y0(20:end)); y1 = ones(n, 1) .* m_end; y1(1:19) = nan;
m_all = mean(y0); x2 = x0; y2 = ones(n, 1) .* m_all;

p(1) = plot(x0, y0); p(2) = plot(x1, y1, 'r--'); p(3) = plot(x2, y2, 'k--');
legend(p, {str_legend_1, 'Среднее 1990 - 2016', 'Среднее 2009 - 2016'})
title(str_title);
saveas(gcf, file_name);