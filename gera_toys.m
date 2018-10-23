%
% Gera problemas "toys"
%

function [] = gera_toys()
    opcoes = {};
    opcoes.sigma = 1e-4;
    opcoes.epsilon = 7.45e-8;
    opcoes.quiet = true;
    opcoes.max_iteracao = 10000;
    opcoes.usar_cauchy = true;
    opcoes.usar_central = true;
    opcoes.usar_linear = true;
    opcoes.usar_pareto = true;

    % Problemas Convexos
    rodadas = 20;
    roda_metodos(rodadas , opcoes , @() fds(50)    , 'FDS');
    roda_metodos(rodadas , opcoes , @() jos1(1000) , 'JOS1');
    roda_metodos(rodadas , opcoes , @lov1          , 'Lov1');
    roda_metodos(rodadas , opcoes , @mop7          , 'MOP7');
    roda_metodos(rodadas , opcoes , @sp1           , 'SP1');

    % Problemas Não Convexos
    roda_metodos(rodadas , opcoes , @ff1          , 'FF1');
    roda_metodos(rodadas , opcoes , @hil1         , 'Hil1');
    roda_metodos(rodadas , opcoes , @() mmr5(100) , 'MMR5');
    roda_metodos(rodadas , opcoes , @mop5         , 'MOP5');
    roda_metodos(rodadas , opcoes , @vu1          , 'VU1');

    % Gráfico
    rodadas = 200;
    gera_pontos_2d(0:0.015:1, @hil1, 'Hil1');
    roda_metodos(rodadas, opcoes, @hil1, 'Hil1', true);

    % Problemas com vários n
    rodadas = 20;
    for n = [200 500 1000 2000]
        opcoes.usar_cauchy = false;
        opcoes.usar_central = true;
        opcoes.usar_linear = false;
        opcoes.usar_pareto = false;
        roda_metodos(rodadas , opcoes , @() fds(n)  , ['FDS_n'  , num2str(n)]);

        opcoes.usar_cauchy = false;
        opcoes.usar_central = false;
        opcoes.usar_linear = true;
        opcoes.usar_pareto = true;
        roda_metodos(rodadas , opcoes , @() mmr5(n) , ['MMR5_n' , num2str(n)]);
    end
end

function [f, Jf, x0, lx, ux] = fds(n)
    f  = @(x) [(1:n) * (x - (1:n)').^4 / n^2;
               exp(ones(1, n) * x / n) + x' * x;
               (1:n) .* (n + 1 - (1:n)) * exp(-x) / (n * (n + 1))];
    Jf = @(x) [(1:n) .* 4 .* (x' - (1:n)).^3 / n^2;
               exp(ones(1, n) * x / n) / n + 2 * x';
               -(1:n) .* (n + 1 - (1:n)) .* exp(-x') / (n * (n + 1))];
    x0 = @() -2 + 4 * rand(n, 1);
    lx = -2 * ones(n, 1);
    ux = 2 * ones(n, 1);
end

function [f, Jf, x0, lx, ux] = jos1(n)
    f  = @(x) [x' * x ; (x - 2)' * (x - 2)] / n;
    Jf = @(x) [x'     ; x' - 2] * 2 / n;
    x0 = @() -1e4 + 2e4 * rand(n, 1);
    lx = [];
    ux = [];
end

function [f, Jf, x0, lx, ux] = lov1()
    f  = @(x) [1.05 * x(1)^2 + 0.98 * x(2)^2;
               0.99 * (x(1) - 3)^2 + 1.03 * (x(2) - 2.5)^2];
    Jf = @(x) [1.05 * 2 * x(1)       , 0.98 * 2 * x(2);
               0.99 * 2 * (x(1) - 3) , 1.03 * 2 * (x(2) - 2.5)];
    x0 = @() -1e2 + 2e2 * rand(2, 1);
    lx = [];
    ux = [];
end

function [f, Jf, x0, lx, ux] = mop7()
    f  = @(x) [(x(1) - 2)^2 / 2 + (x(2) + 1)^2 / 13 + 3;
               (x(1) + x(2) - 3)^2 / 36 + (-x(1) + x(2) + 2)^2 / 8 - 17;
               (x(1) + 2 * x(2) - 1)^2 / 175 + (-x(1) + 2 * x(2))^2 / 17 - 13];
    Jf = @(x) [x(1) - 2                                                      , (x(2) + 1) * 2 / 13;
               (x(1) + x(2) - 3) / 18 - (-x(1) + x(2) + 2) / 4               , (x(1) + x(2) - 3) / 18 + (-x(1) + x(2) + 2) / 4;
               (x(1) + 2 * x(2) - 1) * 2 / 175 - (-x(1) + 2 * x(2)) * 2 / 17 , (x(1) + 2 * x(2) - 1) * 4 / 175 + (-x(1) + 2 * x(2)) * 4 / 17];
    x0 = @() -4e2 + 8e2 * rand(2, 1);
    lx = -4e2 * ones(2, 1);
    ux = 4e2 * ones(2, 1);
end

function [f, Jf, x0, lx, ux] = sp1()
    f  = @(x) [(x(1) - 1)^2 + (x(1) - x(2))^2;
               (x(2) - 3)^2 + (x(1) - x(2))^2];
    Jf = @(x) [(x(1) - 1 + x(1) - x(2)) * 2 , -(x(1) - x(2)) * 2;
               (x(1) - x(2)) * 2            , (x(2) - 3 - (x(1) - x(2))) * 2];
    x0 = @() -1e2 + 2e2 * rand(2, 1);
    lx = [];
    ux = [];
end

function [f, Jf, x0, lx, ux] = ff1()
    f  = @(x) [1 - exp(-(x(1) - 1)^2 - (x(2) + 1)^2);
               1 - exp(-(x(1) + 1)^2 - (x(2) - 1)^2)];
    Jf = @(x) [exp(-(x(1) - 1)^2 - (x(2) + 1)^2) * 2 * (x(1) - 1), exp(-(x(1) - 1)^2 - (x(2) + 1)^2) * 2 * (x(2) + 1);
               exp(-(x(1) + 1)^2 - (x(2) - 1)^2) * 2 * (x(1) + 1), exp(-(x(1) + 1)^2 - (x(2) - 1)^2) * 2 * (x(2) - 1)];
    x0 = @() -1 + 2 * rand(2, 1);
    lx = [];
    ux = [];
end

function [f, Jf, x0, lx, ux] = hil1()
    ac = 45;
    a1 = 40;
    a2 = 25;
    d  = 0.5;
    a  = @(x) (2 * pi / 360) * (ac + [a1 a2] * sin(2 * pi * x));
    da = @(x) (2 * pi / 360) * [a1 a2] .* (2 * pi) .* cos(2 * pi * x');
    b  = @(x) 1 + d * cos(2 * pi * x(1));
    db = @(x) [-2 * pi * d * sin(2 * pi * x(1)), 0];
    f  = @(x) [cos(a(x)) * b(x);
               sin(a(x)) * b(x)];
    Jf = @(x) [-sin(a(x)) * da(x) * b(x) + cos(a(x)) * db(x);
               cos(a(x)) * da(x) * b(x) + sin(a(x)) * db(x)];
    x0 = @() rand(2, 1);
    lx = [];
    ux = [];
end

function [f, Jf, x0, lx, ux] = hil1_grafico()
    [f, Jf, x0] = hil1();
    lx = zeros(2, 1);
    ux = ones(2, 1);
end

function [f, Jf, x0, lx, ux, nome] = mmr5(n)
    f1  = @(x) ((x' * x - ones(1, n) * 10 * cos(2 * pi * x)) / n + 10)^0.25;
    df1 = @(x) 1 / (4 * f1(x)^3) * (2 * x' + 10 * 2 * pi * sin(2 * pi * x')) / n;
    f   = @(x) [f1(x); f1(x - 1.5)];
    Jf  = @(x) [df1(x); df1(x - 1.5)];
    x0  = @() -5 + 10 * rand(n, 1);
    lx  = -5 * ones(n, 1);
    ux  = 5 * ones(n, 1);
end

function [f, Jf, x0, lx, ux] = mop5()
    f  = @(x) [x' * x / 2 + sin(x' * x);
               (3 * x(1) - 2 * x(2) + 4)^2 / 8 + (x(1) - x(2) + 1)^2 / 27 + 15;
               1 / (x' * x + 1) - 1.1 * exp(-x' * x)];
    Jf = @(x) [(1 + 2 * cos(x' * x)) * x';
               (3 * x(1) - 2 * x(2) + 4) * 3 / 4 + (x(1) - x(2) + 1) * 2 / 27, -(3 * x(1) - 2 * x(2) + 4) / 2 - (x(1) - x(2) + 1) * 2 / 27;
               (-1 / (x' * x + 1)^2 + 1.1 * exp(-x' * x)) * 2 * x'];
    x0 = @() -1 + 2 * rand(2, 1);
    lx = -30 * ones(2, 1);
    ux = 30 * ones(2, 1);
end

function [f, Jf, x0, lx, ux] = vu1()
    f  = @(x) [1 / (x' * x + 1);
               x(1)^2 + 3 * x(2)^2 + 1];
    Jf = @(x) [-1 / (x' * x + 1)^2 * 2 * x';
               2 * x(1), 6 * x(2)];
    x0 = @() -3 + 6 * rand(2, 1);
    lx = -3 * ones(2, 1);
    ux = 3 * ones(2, 1);
end

function [] = gera_pontos_2d(malha, problema, nome)
    saida = ['gerado/toy_' nome '_pontos.dat'];
    if ~exist(saida, 'file')
        mostra([cor(37, 1) 'Gerando ' saida cor()]);
    else
        mostra([cor(37, 2) 'Pulando ' saida cor()]);
        return
    end
    f = problema();
    file = fopen(saida, 'w');
    for x = malha
        for y = malha
            fprintf(file, [num2str(f([x; y])') '\n']);
        end
    end
    fclose(file);
end
function [] = roda_metodos(rodadas, opcoes, problema, nome, grafico)
    if nargin < 5
        grafico = false;
    end

    cauchy = {};
    cauchy.saida = ['gerado/toy_' nome '_grafico_cauchy.dat'];
    cauchy.sucesso = 0;
    cauchy.evalf = [];
    cauchy.tempo = [];
    cauchy.iteracoes = [];
    cauchy.xs = [];
    cauchy.x0s = [];

    central = {};
    central.saida = ['gerado/toy_' nome '_grafico_central.dat'];
    central.sucesso = 0;
    central.evalf = [];
    central.tempo = [];
    central.iteracoes = [];
    central.xs = [];
    central.x0s = [];

    linear = {};
    linear.saida = ['gerado/toy_' nome '_grafico_linear.dat'];
    linear.sucesso = 0;
    linear.evalf = [];
    linear.tempo = [];
    linear.iteracoes = [];
    linear.xs = [];
    linear.x0s = [];

    pareto = {};
    pareto.saida = ['gerado/toy_' nome '_grafico_pareto.dat'];
    pareto.sucesso = 0;
    pareto.evalf = [];
    pareto.tempo = [];
    pareto.iteracoes = [];
    pareto.xs = [];
    pareto.x0s = [];

    if grafico
        saida = pareto.saida;
    else
        saida = ['gerado/toy_' nome '.txt'];
    end
    if ~exist(saida, 'file')
        mostra([cor(37, 1) 'Gerando ' saida cor()]);
    else
        mostra([cor(37, 2) 'Pulando ' saida cor()]);
        return
    end

    [f, Jf, x0, lx, ux] = problema();

    for r = 1:rodadas
        mostra(['Rodada ' num2str(r)]);
        x0r = x0();

        if opcoes.usar_cauchy
            mostra('Cauchy');
            tic;
            opcoes.criterio_parada = 2;
            [xsol, fsol, k, flag, dados] = metodo_descida(f, @direcao_cauchy, Jf, x0r, lx, ux, [], [], [], [], opcoes);
            if flag == 1
                if grafico
                    cauchy.xs = [cauchy.xs xsol];
                    cauchy.x0s = [cauchy.x0s x0r];
                else
                    cauchy.sucesso = cauchy.sucesso + 1;
                    cauchy.evalf = [cauchy.evalf, dados.evalf];
                    cauchy.tempo = [cauchy.tempo, dados.tempo];
                    cauchy.iteracoes = [cauchy.iteracoes, k];
                end
            end
            toc;
        end

        if opcoes.usar_central
            mostra('Central');
            tic;
            opcoes.criterio_parada = 2;
            [xsol, fsol, k, flag, dados] = metodo_descida(f, @direcao_central, Jf, x0r, lx, ux, [], [], [], [], opcoes);
            if flag == 1
                if grafico
                    central.xs = [central.xs xsol];
                    central.x0s = [central.x0s x0r];
                else
                    central.sucesso = central.sucesso + 1;
                    central.evalf = [central.evalf, dados.evalf];
                    central.tempo = [central.tempo, dados.tempo];
                    central.iteracoes = [central.iteracoes, k];
                end
            end
            toc;
        end

        if opcoes.usar_linear
            mostra('Linear');
            tic;
            opcoes.criterio_parada = 1;
            [xsol, fsol, k, flag, dados] = metodo_descida(f, @direcao_linear, Jf, x0r, lx, ux, [], [], [], [], opcoes);
            if flag == 1
                if grafico
                    linear.xs = [linear.xs xsol];
                    linear.x0s = [linear.x0s x0r];
                else
                    linear.sucesso = linear.sucesso + 1;
                    linear.evalf = [linear.evalf, dados.evalf];
                    linear.tempo = [linear.tempo, dados.tempo];
                    linear.iteracoes = [linear.iteracoes, k];
                end
            end
            toc;
        end

        if opcoes.usar_pareto
            mostra('Pareto');
            tic;
            opcoes.criterio_parada = 1;
            [xsol, fsol, k, flag, dados] = metodo_descida(f, @direcao_pareto, Jf, x0r, lx, ux, [], [], [], [], opcoes);
            if flag == 1
                if grafico
                    pareto.xs = [pareto.xs xsol];
                    pareto.x0s = [pareto.x0s x0r];
                else
                    pareto.sucesso = pareto.sucesso + 1;
                    pareto.evalf = [pareto.evalf, dados.evalf];
                    pareto.tempo = [pareto.tempo, dados.tempo];
                    pareto.iteracoes = [pareto.iteracoes, k];
                end
            end
            toc;
        end
    end

    if grafico
        imprime_saida(f, cauchy.x0s, cauchy.xs, cauchy.saida);
        imprime_saida(f, central.x0s, central.xs, central.saida);
        imprime_saida(f, linear.x0s, linear.xs, linear.saida);
        imprime_saida(f, pareto.x0s, pareto.xs, pareto.saida);
    else
        sucessos = [cauchy.sucesso central.sucesso linear.sucesso pareto.sucesso];
        evalfs = [mediana(cauchy.evalf) mediana(central.evalf) mediana(linear.evalf) mediana(pareto.evalf)];
        tempos = [mediana(cauchy.tempo) mediana(central.tempo) mediana(linear.tempo) mediana(pareto.tempo)];
        iteracoes = [mediana(cauchy.iteracoes) mediana(central.iteracoes) mediana(linear.iteracoes) mediana(pareto.iteracoes)];

        file = fopen(saida, 'w');
        mostraSaida(file,  '          cauchy central  linear  pareto');
        mostraSaida(file, ['Suc.(%):' formata(sucessos * 100 / rodadas)]);
        mostraSaida(file, ['Iter.:  ' formata(iteracoes)]);
        mostraSaida(file, ['Eval.f: ' formata(evalfs)]);
        mostraSaida(file, ['Tem.(s):' formata(tempos)]);
        mostra;
        fclose(file);
    end
end

function [] = imprime_saida(f, x0s, xs, saida)
    file = fopen(saida, 'w');
    for k = 1:2:size(xs, 2)
        fprintf(file, [num2str(f(x0s(:, k))') ' a\n']);
        fprintf(file, [num2str(f(xs(:, k))') ' b\n']);
        fprintf(file, '\n');
    end
    fclose(file);
end

function [] = mostraSaida(file, msg)
    mostra(msg);
    fprintf(file, '%s\n', msg);
end

function m = mediana(arr)
    if length(arr) > 0
        m = median(arr);
    else
        m = -1;
    end
end

function s = formata(arr)
    s = '';
    for k = 1:length(arr)
        if arr(k) == -1
            s = [s sprintf(' %7s', '-')];
        else
            s = [s sprintf(' %7.1f', arr(k))];
        end
    end
end
