%
%
%
% Saída:
% ['gerado/portifolio_' portifolio '.txt'] : resultados com formatação de cor e a cada rodada
%

function [] = gera_portifolio()
    opcoes = {};
    mostra;
    mostra([cor(36, 1) 'Gerando exemplos de portifolio' cor()]);
    metodos_comparados = [ % precisa ter exatamente 7 caracteres cada
        '  CZeSD';
        ' KP_SSD';
        '  L_SSD';
        'LR_ASSD';
        'MeanVar';
        'RMZ_SSD';
    ];
    opcoes.epsilon = 1e-7;
    opcoes.sigma = 1e-4;
    opcoes.quiet = true;
    portifolio('DowJones', metodos_comparados, opcoes);       % 28 assets, 110 rebalancing
    portifolio('FF49Industries', metodos_comparados, opcoes); % 49 assets, 190 rebalancing
    portifolio('NASDAQ100', metodos_comparados, opcoes);      % 82 assets,  46 rebalancing
    portifolio('FTSE100', metodos_comparados, opcoes);        % 83 assets,  56 rebalancing

    opcoes.epsilon = 1e-6;
    portifolio('SP500', metodos_comparados, opcoes);          %  442 assets, 46 rebalancing
    portifolio('NASDAQComp', metodos_comparados, opcoes);     % 1203 assets, 53 rebalancing
end

% gera arquivo de saída, onde cada linha é o retorno total obtido no outsample. Cada coluna é de um método diferente (pareto, linear, etc)
function [] = portifolio(nome, metodos_comparados, opcoes)
    saida = ['gerado/portifolio_' nome '.txt'];
    if exist(saida, 'file')
        mostra([cor(37, 2) 'Pulando ' saida cor()]);
        return
    else
        mostra([cor(37, 1) 'Gerando ' saida cor()]);
    end

    % inicializa
    load(['entrada/' nome '.mat']);

    [samples, n] = size(Assets_Returns);

    A = ones(1, n);
    a = 1;
    lx = zeros(n, 1);
    x0 = ones(n, 1) / n;

    undef = -9999;
    metodosDescida = 4;

    qtd_comparados = size(metodos_comparados, 1);
    comparados = [];
    for c = 1:qtd_comparados
        variavel = ['OptPortfolios_' strtrim(metodos_comparados(c, :)) '_' nome];
        load(['entrada/' variavel '.txt']);
        comparados(:, :, c) = eval(variavel);
    end

    resultadosAcumulados = zeros(1, qtd_comparados + metodosDescida);
    medicoesAcumuladas = zeros(1, 2 * metodosDescida);
    sucessosAcumulados = zeros(1, metodosDescida);

    file = fopen(saida, 'w');

    rodadas = ceil((samples - 52) / 12);

    cabecalho  = [cor(37, 1) '  cauchy central  linear  pareto'];
    for c = 1:qtd_comparados
         cabecalho  = [cabecalho  ' ' metodos_comparados(c, :)];
    end
    cabecalho  = [cabecalho  '       cauchy      central       linear       pareto' cor()];

    fprintf(file, [cor(37, 1) nome cor() '\n']);
    fprintf(file, [cabecalho '\n']);

    % loop onde cada rodada equivale a obter o portifolio que minimiza o risco e maximiza o retorno das últimas 52 semanas, e aplica nas próximas 12.
    % O loop seguinte começa 12 semanas depois, portanto parte das 52 semanas ainda são usadas.
    for t = 1:rodadas
        % 52 semanas de dados para usar para decidir os assets
        ti = (t - 1) * 12 + 1;
        tf = ti + 51;
        insample = Assets_Returns(ti:tf, :);

        % 12 semanas de dados para testar se os assets decididos foram bons
        ti = tf + 1;
        tf = min(samples, ti + 11);
        outsample = Assets_Returns(ti:tf, :);

        % gera a matriz de covariância
        C = cov(insample);
        if n > 500
            C = sparse(C);
        end

        % minimiza a variância do portifólio e o retorno médio
        f = @(x) [-ones(1, 52) * insample * x; 0.5 * x' * C * x];
        Jf = @(x) [-(ones(1, 52) * insample)', C * x]';

        medicoes = [];
        resultados = [];
        sucessos = [];

        [xsol, iteracoes, tempo] = roda_portifolio(undef, opcoes, @(op) metodo_descida(f, @direcao_cauchy, Jf, x0, lx, [], A, a, [], [], op));
        resultados = [resultados retornoTotal(undef, outsample, xsol)];
        medicoes = [medicoes iteracoes tempo];
        sucessos = [sucessos ~isempty(xsol)];

        [xsol, iteracoes, tempo] = roda_portifolio(undef, opcoes, @(op) metodo_descida(f, @direcao_central, Jf, x0, lx, [], A, a, [], [], op));
        resultados = [resultados retornoTotal(undef, outsample, xsol)];
        medicoes = [medicoes iteracoes tempo];
        sucessos = [sucessos ~isempty(xsol)];

        [xsol, iteracoes, tempo] = roda_portifolio(undef, opcoes, @(op) metodo_descida(f, @direcao_linear, Jf, x0, lx, [], A, a, [], [], op));
        resultados = [resultados retornoTotal(undef, outsample, xsol)];
        medicoes = [medicoes iteracoes tempo];
        sucessos = [sucessos ~isempty(xsol)];

        [xsol, iteracoes, tempo] = roda_portifolio(undef, opcoes, @(op) metodo_descida(f, @direcao_pareto, Jf, x0, lx, [], A, a, [], [], op));
        resultados = [resultados retornoTotal(undef, outsample, xsol)];
        medicoes = [medicoes iteracoes tempo];
        sucessos = [sucessos ~isempty(xsol)];

        for c = 1:qtd_comparados
            resultados = [resultados retornoTotal(undef, outsample, comparados(:, t, c))];
        end

        imprimeResultados(undef, resultados, medicoes, file);

        resultadosAcumulados(resultados ~= undef) = resultadosAcumulados(resultados ~= undef) + resultados(resultados ~= undef);
        medicoesAcumuladas(medicoes ~= undef) = medicoesAcumuladas(medicoes ~= undef) + medicoes(medicoes ~= undef);
        sucessosAcumulados = sucessosAcumulados + sucessos;
    end
    fprintf(file, [cor(37, 1) 'Média:' cor() '\n']);
    fprintf(file, [cabecalho '\n']);
    imprimeResultados(undef, resultadosAcumulados, medicoesAcumuladas, file, sucessosAcumulados, rodadas);
    fclose(file);
end

function [] = imprimeResultados(undef, resultados, medicoes, file, sucessos, rodadas)
    metodosDescida = length(medicoes) / 2;
    formato = '';

    if nargin >= 5
        for m = 1:metodosDescida
            sucesso = sucessos(m);
            if sucesso == 0
                resultados(m) = undef;
                medicoes(2 * m - 1) = undef;
                medicoes(2 * m) = undef;
            else
                resultados(m) = resultados(m) / sucesso;
                medicoes(2 * m - 1) = medicoes(2 * m - 1) / sucesso;
                medicoes(2 * m) = medicoes(2 * m) / sucesso;
            end
        end
        resultados((metodosDescida + 1):end) = resultados((metodosDescida + 1):end) / rodadas;
    end

    resultadosFiltrados = resultados(resultados ~= undef);
    resultadosDescida = resultados(1:metodosDescida);
    resultadosDescidaFiltrados = resultadosDescida(resultadosDescida ~= undef);
    resultadosOutros = resultados((metodosDescida + 1) : end);
    resultadosOutrosFiltrados = resultadosOutros(resultadosOutros ~= undef);
    for r = 1:length(resultados)
        if resultados(r) == max(resultadosFiltrados)
            corAdotada = cor(36, 1);
        elseif resultados(r) == min(resultadosFiltrados)
            corAdotada = cor(37, 2);
        elseif r <= metodosDescida
            if resultados(r) == max(resultadosDescidaFiltrados)
                corAdotada = cor(32, 1);
            else
                corAdotada = cor(32);
            end
        else
            if resultados(r) == max(resultadosOutrosFiltrados)
                corAdotada = cor(33, 1);
            else
                corAdotada = cor(33);
            end
        end
        if resultados(r) == undef
            formato = [formato ' ' corAdotada '   -   ' cor()];
        else
            formato = [formato ' ' corAdotada '%7.3f' cor()];
        end
    end

    formato = [formato  ' '];

    for m = 1:metodosDescida
        if medicoes(2 * m) == undef
            formato = [formato cor(31) '  -  ' cor(34, 1) '   -    '];
        else
            formato = [formato cor(31) '%4.0f ' cor(34, 1) '%6.1fs '];
        end
    end
    formato = [formato cor() '\n'];

    medicoessFiltradas = medicoes(medicoes ~= undef);
    fprintf(file, formato, [resultadosFiltrados medicoessFiltradas]);
    fflush(file);
end

function resultado = retornoTotal(undef, outsample, x)
    if isempty(x)
        resultado = undef;
    else
        resultado = sum(outsample * x);
    end
end

function [xsol, k, tempo] = roda_portifolio(undef, opcoes, exemplo)
    tic;
    [xsol, fsol, k, flag] = exemplo(opcoes);
    tempo = toc;
    if flag ~= 1 && flag ~= 2
        xsol = [];
        tempo = undef;
        k = undef;
    end
end


