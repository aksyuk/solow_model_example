% 23 Jan 2020
% Светлана Суязова (Аксюк)
% Модель Солоу на данных Всемирного банка по экономике США c 1990 по 2016


% Загрузить статистику World Bank по стране _______________________________
%   показатели:
%   * NE.GDI.FTOT.KD      Gross fixed capital formation (constant 2000 US$)
%   * NY.GDP.MKTP.KD      GDP (constant 2000 US$)
%   * SL.TLF.TOTL.IN      Labor force, total
%   * SP.POP.TOTL         Population, total
%   * GC.NFN.TOTL.GD.ZS   Net investment in nonfinancial assets (% of GDP)

% скрипт, сгенерированный процедурой импорта ------------------------------
% Setup the Import Options and import the data
opts = delimitedTextImportOptions("NumVariables", 9);

% Specify range and delimiter
opts.DataLines = [2, 29];
opts.Delimiter = ";";

% Specify column names and types
opts.VariableNames = ["iso2c", "country", "year", "NE_GDI_FTOT_KD", ...
    "NY_GDP_MKTP_KD", "GC_NFN_TOTL_GD_ZS", "SL_TLF_TOTL_IN", ...
    "SP_POP_TOTL", "NY_GNS_ICTR_ZS"];
opts.VariableTypes = ["categorical", "categorical", "double", "double", ... 
    "double", "double", "double", "double", "double"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Specify variable properties
opts = setvaropts(opts, ["iso2c", "country"], "EmptyFieldRule", "auto");
opts = setvaropts(opts, ["year", "NE_GDI_FTOT_KD", "NY_GDP_MKTP_KD", ...
    "GC_NFN_TOTL_GD_ZS", "SL_TLF_TOTL_IN", "SP_POP_TOTL", ...
    "NY_GNS_ICTR_ZS"], "DecimalSeparator", ",");

% Import the data
dt = readtable("~/Documents/MATLAB/solow_data_US.csv", opts);

% Clear temporary variables
clear opts


% Загрузить статистику амортизации капитала по США ________________________
% источник: https://apps.bea.gov/iTable/iTable.cfm?reqid=19&step=2
%  Table 1.7.5. Relation of Gross Domestic Product, 
%  Gross National Product, Net National Product, National Income, 
%  and Personal Income
% Consumption of fixed capital, Billions of dollars (1990 - 2016, annual)

% скрипт, сгенерированный процедурой импорта ------------------------------
% Setup the Import Options and import the data
opts = delimitedTextImportOptions("NumVariables", 29);

% Specify range and delimiter
opts.DataLines = [5, 5; 7, 7; 11, 11];
opts.Delimiter = ",";

% Specify column names and types
opts.VariableNames = ["Var1", "Var2", "VarName3", "VarName4", ...
    "VarName5", "VarName6", "VarName7", "VarName8", "VarName9", ...
    "VarName10", "VarName11", "VarName12", "VarName13", "VarName14", ...
    "VarName15", "VarName16", "VarName17", "VarName18", "VarName19", ...
    "VarName20", "VarName21", "VarName22", "VarName23", "VarName24", ...
    "VarName25", "VarName26", "VarName27", "VarName28", "VarName29"];
opts.SelectedVariableNames = ["VarName3", "VarName4", "VarName5", ...
    "VarName6", "VarName7", "VarName8", "VarName9", "VarName10", ...
    "VarName11", "VarName12", "VarName13", "VarName14", "VarName15", ...
    "VarName16", "VarName17", "VarName18", "VarName19", "VarName20", ...
    "VarName21", "VarName22", "VarName23", "VarName24", "VarName25", ...
    "VarName26", "VarName27", "VarName28", "VarName29"];
opts.VariableTypes = ["string", "string", "double", "double", "double", ...
    "double", "double", "double", "double", "double", "double", ...
    "double", "double", "double", "double", "double", "double", ...
    "double", "double", "double", "double", "double", "double", ...
    "double", "double", "double", "double", "double", "double"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Specify variable properties
opts = setvaropts(opts, ["Var1", "Var2"], "WhitespaceRule", "preserve");
opts = setvaropts(opts, ["Var1", "Var2"], "EmptyFieldRule", "auto");

% Import the data
depr = readtable("/home/light/Documents/MATLAB/bea_table_1-7-5.csv", opts);

% Clear temporary variables
clear opts


% Оцениваем параметры модели ______________________________________________
% переменная I: инвестиции, млн. долл. США
I = dt.NY_GDP_MKTP_KD .* dt.GC_NFN_TOTL_GD_ZS / 100 / 10^6;
% переменная Y: ВВП, млн. долл. США
Y = dt.NY_GDP_MKTP_KD ./ 10^6;
% переменная K: стоимость ОПФ, млн. долл. США
K = dt.NE_GDI_FTOT_KD ./ 10^6;
% переменная L: численность рабочей силы, млн. чел.
L = dt.SL_TLF_TOTL_IN ./ 10^6;

% функция Кобба-Дугласа (МНК на натуральных логарифмах)
X = [ones(length(K),1) log(K), log(L)];
b = X \ log(Y);
Y_fit = X * b;

% ищем параметры нелинейной целевой функции нескольких переменных
x0 = [0.0001, 0.0001, 0.0001];
n = length(Y);
count = 1;

res = array2table(zeros(n-1, 5), 'VariableNames', {'w', 'A', 'alpha_K', ...
    'alpha_L', 'err'});
for i = 2:n
    T = @(x) CB_target(K, L, Y, x, i);
    [x2, fval2] = fmincon(T, x0, [-1, 0, 0], 0, [0, 1, 1], 1, ...
        [0, 0, 0], [inf, 1, 1]);
    res{count, 1} = i; res{count, 2:4} = x2; res{count, 5} = fval2;
    count = count + 1;
end
indx = find(res.err(2:end) - res.err(1:end-1) == ...
    min(res.err(2:end) - res.err(1:end-1)));
x2 = res{indx + 1, 2:4};
Y_fit = x2(1) .* K .^ x2(2) .* L .^ x2(3);

% график прогноз-реализация
plot_obs_vs_pred(Y, Y_fit)

x2 = res{end, 2:4};
Y_fit = x2(1) .* K .^ x2(2) .* L .^ x2(3);

hold on
scatter(Y, Y_fit)
legend('Модель с окном 25 наблюдений', 'Идеальный прогноз', ...
    'Модель по всем 27 наблюдениям')
saveas(gcf, "~/Documents/MATLAB/plots/plot_mod_vs_fit.png");


% параметры модели --------------------------------------------------------
%  * фактор шкалы
% A = exp(b(1));      % по МНК, не функция Кобба-Дугласа
A = res{end, "A"};

%  * эластичности выпуска по капиталу и труду
% alpha_K = b(2);     % по МНК, не функция Кобба-Дугласа
% alpha_L = b(3);     % по МНК, не функция Кобба-Дугласа
alpha_K = res{end, "alpha_K"}; 
alpha_L = res{end, "alpha_L"};

%  * норма инвестирования
x0 = dt.year; y0 = dt.NY_GNS_ICTR_ZS;
[mean_all, mean_end] = plot_means(x0, y0, 'Норма инвестирования, %', ...
    'Динамика нормы инвестирования в США', ...
    "~/Documents/MATLAB/plots/plot_means-01.png");
%  считаем по данным после кризиса
rho = mean_end / 100;

%  * норма амортизации
x0 = depr{1, :}; y0 = depr{3, :} ./ depr{2, :} * 100;
[mean_all, mean_end] = plot_means(x0, y0, 'Норма амортизации, %', ...
    'Динамика нормы амортизации в США', ...
    "~/Documents/MATLAB/plots/plot_means-02.png");
%  считаем по данным после кризиса
mu = mean_end / 100; 

%  * темп прироста рабочей силы
x0 = dt.year(2:end);
y0 = dt.SL_TLF_TOTL_IN(2:end) ./ dt.SL_TLF_TOTL_IN(1:end-1) - 1; 
[mean_all, mean_end] = plot_means(x0, y0, 'Темп прироста рабочей силы', ...
    'Динамика темпа прироста рабочей силы в США', ...
    "~/Documents/MATLAB/plots/plot_means-03.png");
%  считаем по данным после кризиса
nu = mean_end;  


% Переходный режим в модели Солоу _________________________________________
[vt, v0] = solow_time_table(K(end), L(end), A, alpha_K, alpha_L, rho, ...
    mu, nu, 300, "~/Documents/MATLAB/plots/plot_trace-01.png");

% график k vs y в точке равновесия
clf
n = length(vt.t);
m = 50;
x0 = linspace(0, v0{1, 'k'} * 2, m);
y0 = A .* x0 .^ alpha_K; y1 = y0 * (1-rho); y2 = x0 .* (mu+nu)*(1-rho)/rho;
plot(x0, y0, 'k', x0, y1, 'r', x0, y2, 'b', ...
    v0{1, 'k'} .* ones(m), linspace(0, v0{1, 'y'}, m), 'k--', ...
    linspace(0, v0{1, 'k'}, m), v0{1, 'y'} .* ones(m), 'k--', ...
    linspace(0, v0{1, 'k'}, m), v0{1, 'y'} * (1-rho) .* ones(m), 'k--')
xlabel('k'); ylabel('y'); title('Модель Солоу, равновесие');
legend({'\it{f(k)}', '\it{(1-\rho)f(k)}', '\it{k(\mu+\nu)(1-\rho)/\rho}'});
text([0, 0, v0{1, 'k'}] + 10^5, ...
    [v0{1, 'y'}, v0{1, 'y'} * (1-rho), 0] + 10^6/2, ...
    {'\it{f(k^*)}', '\it{(1-\rho)f(k^*)}', '\it{k^*}'});
saveas(gcf, "~/Documents/MATLAB/plots/plot_eq-01.png");

% % сохраняем данные
% %  * для переменных
% tbl = array2table([dt.year, K, L, Y, I], 'VariableNames', ...
%     {'год', 'K', 'L', 'Y', 'I'});
% writetable(tbl, "~/Documents/MATLAB/solow_vars_matlab.csv", ...
%     'Delimiter', ';');
% %  * для параметров: норма амортизации
% depr.Properties.VariableNames = string(depr{1, :}); depr(1, :) = [];
% writetable(depr, "~/Documents/MATLAB/solow_depr_matlab.csv", ...
%     'Delimiter', ';');
%  * для параметров: норма инвестирования и темп прироста рабочей силы
% mu_stats = [NaN; dt.SL_TLF_TOTL_IN(2:end) ./ dt.SL_TLF_TOTL_IN(1:end-1)-1];
% tbl = array2table([dt.year, dt.NY_GNS_ICTR_ZS, mu_stats], 'VariableNames', ...
%     {'год', 'rho', 'nu'});
% writetable(tbl, "~/Documents/MATLAB/solow_param_matlab.csv", ...
%     'Delimiter', ';');
