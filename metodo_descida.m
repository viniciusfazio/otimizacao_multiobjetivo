%
% Método de Descida Multiobjetivo com backtracking
%
% Entrada:
% f                           : função R^n -> R^m a ser minimizada
% vf                          : função que retorna a direção de descida. function [v, msg_v, flag] = vf(x, k, J, Jf, lx, ux, A, a, B, b, opcoes)
%                             : onde x      = ponto iterado do R^n
%                             :      k      = número da iteração
%                             :      J      = matriz jacobiana do ponto x, ou seja, Jf(x)
%                             :      Jf     = parâmetro Jf descrito abaixo
%                             :      lx     = valor mínimo para x + v. [] significa sem limite mínimo.
%                             :      ux     = valor máximo para x + v. [] significa sem limite máximo.
%                             :      A      = restrição linear A(x + v) = a. Deve ser compatível com a. [] significa sem restrição.
%                             :      a      = restrição linear A(x + v) = a. Deve ser compatível com A. [] significa sem restrição.
%                             :      B      = restrição linear B(x + v) <= b. Deve ser compatível com b. [] significa sem restrição.
%                             :      b      = restrição linear B(x + v) <= b. Deve ser compatível com B. [] significa sem restrição.
%                             :      opcoes = parâmetro opcoes descrito abaixo
%                             :      v      = vetor R^n direção de descida. Caso v não seja de descida, o algoritmo para.
%                             :      msg_v  = string contendo alguma informação para mostrar na tela a cada iteração
%                             :      flag   = 1: sucesso
%                             :               2: erro
% Jf                          : função que retorna a matriz jacobiana de um ponto qualquer. function J = Jf(x)
%                             : onde x = ponto do R^n
%                             :      J = matriz m x n jacobiana do ponto x
% x0                          : ponto inicial do R^n
% lx                          : (padrão []) valor mínimo para x. [] significa sem limite mínimo.
% ux                          : (padrão []) valor máximo para x. [] significa sem limite máximo.
% A                           : (padrão []) restrição linear Ax = a. Deve ser compatível com a. [] significa sem restrição.
% a                           : (padrão []) restrição linear Ax = a. Deve ser compatível com A. [] significa sem restrição.
% B                           : (padrão []) restrição linear Bx <= b. Deve ser compatível com b. [] significa sem restrição.
% b                           : (padrão []) restrição linear Bx <= b. Deve ser compatível com B. [] significa sem restrição.
% opcoes                      : (padrão {}) estrutura contendo as demais opções
% opcoes.callback             : (padrão []) função chamada no início de cada iteração. function [] = callback(x, k)
%                             : onde x = ponto iterado do R^n
%                             :      k = número da iteração
% opcoes.criterio_parada      : (padrão 1) 1: max Jv >= -epsilon
%                             :            2: max Jv + 0.5|v|^2 >= -epsilon
% opcoes.mostra_cada_iteracao : (padrão false) imprime um resumo de cada iteração enquanto roda
% opcoes.max_iteracao         : (padrão 10000) quantidade máxima de iterações
% opcoes.max_tempo            : (padrão 9999999) quantidade máxima de segundos que o método pode demorar
% opcoes.sigma                : (padrão 0.5) constante do Método de Descida, no intervalo (0, 1)
% opcoes.epsilon              : (padrão 1e-8) epsilon usado no critério de parada
% opcoes.quiet                : (padrão false) se true, não imprime nenhuma mensagem
%
% Saída:
% xsol           : ponto solução
% fsol           : f(xsol)
% k              : número de iterações
% flag           : 1 = sucesso
%                : 2 = máximo de iterações ou tempo máximo atingido
%                : 3 = erro nos parâmetros
%                : 4 = erro no backtracking
%                : 5 = erro na direção de descida
% dados          : estrutura contendo todos dados das iterações
% dados.evalf    : quantidade de vezes que a função f foi chamada no total
% dados.fs       : f aplicado cada ponto iterados
% dados.sigmaJvs : todos sigma * J * v usado para obter o t
% dados.tempo    : tempo decorrido no total
% dados.ts       : todos escalares t obtidos
% dados.vs       : todas direções obtidas
% dados.xs       : todos os pontos iterados
%

function [xsol, fsol, k, flag, dados] = metodo_descida(f, vf, Jf, x0, lx, ux, A, a, B, b, opcoes)
    % Constantes
    gerar_dados = nargout >= 5;

    % Inicializa saída
    xsol = 0;
    fsol = 0;
    k    = 0;
    flag = 3;
    if gerar_dados
        dados          = {};
        dados.evalf    = 0;
        dados.fs       = [];
        dados.sigmaJvs = [];
        dados.tempo    = 0;
        dados.ts       = [];
        dados.vs       = [];
        dados.xs       = [];
    end

    % Verifica se os parâmetros foram passados
    if nargin < 4
        mostra('Quantidade de parâmetros inválido. É necessário informar f, vf, Jf, x0');
        return
    end
    if nargin < 5
        lx = [];
    end
    if nargin < 6
        ux = [];
    end
    if nargin < 7
        A = [];
    end
    if nargin < 8
        a = [];
    end
    if nargin < 9
        B = [];
    end
    if nargin < 10
        b = [];
    end
    if nargin < 11
        opcoes = {};
    end
    if ~isfield(opcoes, 'callback')
        opcoes.callback = [];
    end
    if ~isfield(opcoes, 'mostra_cada_iteracao')
        opcoes.mostra_cada_iteracao = false;
    end
    if ~isfield(opcoes, 'max_iteracao')
        opcoes.max_iteracao = 10000;
    end
    if ~isfield(opcoes, 'max_tempo')
        opcoes.max_tempo = 9999999;
    end
    if ~isfield(opcoes, 'criterio_parada')
        opcoes.criterio_parada = 1;
    end
    if ~isfield(opcoes, 'sigma')
        opcoes.sigma = 0.5;
    end
    if ~isfield(opcoes, 'epsilon')
        opcoes.epsilon = 1e-8;
    end
    if ~isfield(opcoes, 'quiet')
        opcoes.quiet = false;
    end
    tem_callback = isa(opcoes.callback, 'function_handle');

    % Verifica se os parâmetros são válidos
    if ~isa(f, 'function_handle')
        mostra('f deve ser uma função.');
        return
    end

    if ~isa(vf, 'function_handle')
        mostra('vf deve ser uma função.');
        return
    end

    if ~isa(Jf, 'function_handle')
        mostra('Jf deve ser uma função.');
        return
    end

    if ~isnumeric(opcoes.sigma) || opcoes.sigma <= 0 || opcoes.sigma >= 1
        mostra('opcoes.sigma deve ser um número entre 0 e 1, sem incluir 0 e 1.');
        return
    end

    if ~isnumeric(opcoes.max_iteracao) || opcoes.max_iteracao < 1
        mostra('opcoes.max_iteracao deve ser um número maior ou igual a 1.');
        return
    end

    if ~isnumeric(opcoes.max_tempo) || opcoes.max_tempo < 1
        mostra('opcoes.max_tempo deve ser um número maior ou igual a 1.');
        return
    end

    if ~isnumeric(opcoes.criterio_parada) || ~any(opcoes.criterio_parada == [1, 2])
        mostra('opcoes.criterio_parada deve ser 1 ou 2.');
        return
    end

    if ~ismatrix(x0) || size(x0, 2) ~= 1
        mostra('x0 deve ser um vetor coluna de n linhas');
        return
    end
    nx = size(x0, 1);

    fx0 = f(x0);
    if ~ismatrix(fx0) || size(fx0, 2) ~= 1
        mostra('f deve ser uma função vetorial, que retorna um vetor coluna de m linhas');
        return
    end
    mf = size(fx0, 1);

    Jfx0 = Jf(x0);
    if ~ismatrix(Jfx0) || any(size(Jfx0) ~= [mf, nx])
        mostra('Jf deve retornar uma matriz de m linhas e n colunas, onde m é a quantidade de linhas do vetor coluna f(x0) e n é a quantidade de linhas do vetor coluna x0.');
        return
    end
    [m, n] = size(Jfx0);

    if ~ismatrix(lx) || ~isempty(lx) && any(size(lx) ~= [n, 1])
        mostra('lx deve ser um vetor coluna com n linhas, ou ser vazio.');
        return
    end

    if ~ismatrix(ux) || ~isempty(ux) && any(size(ux) ~= [n, 1])
        mostra('ux deve ser um vetor coluna com n linhas, ou ser vazio.');
        return
    end

    if ~ismatrix(A) || ~isempty(A) && size(A, 2) ~= n
        mostra('A matriz A deve ter n colunas, ou ser vazio.');
        return
    end

    if ~ismatrix(a) || ~isempty(A) && any(size(a) ~= [size(A, 1), 1])
        mostra('O vetor a deve ser um vetor coluna com a mesma quantidade de linhas que a matriz A, ou ser vazio.');
        return
    end

    if ~ismatrix(B) || ~isempty(B) && size(B, 2) ~= n
        mostra('A matriz B deve ter n colunas, ou ser vazio.');
        return
    end

    if ~ismatrix(b) || ~isempty(B) && any(size(b) ~= [size(B, 1), 1])
        mostra('O vetor b deve ser um vetor coluna com a mesma quantidade de linhas que a matriz B, ou ser vazio.');
        return
    end

    % início das iterações
    xsol = x0;
    fsol = f(xsol);
    if gerar_dados
        dados.fs = [dados.fs, fsol];
        dados.xs = [dados.xs, xsol];
        dados.evalf = 0;
    end
    tic_id = tic;
    for k = 0:(opcoes.max_iteracao - 1)
        if toc(tic_id) >= opcoes.max_tempo
            k = k - 1;
            break
        end
        if tem_callback
            opcoes.callback(xsol, k);
        end
        tic_iter_id = tic;
        J = Jf(xsol);
        [v, msg_v, flag_direcao] = vf(xsol, k, J, Jf, lx, ux, A, a, B, b, opcoes);
        if any(size(v) ~= [n, 1])
            if gerar_dados
                dados.tempo = toc(tic_id);
            end
            mostra(['Direção de descida com dimensões inválidas: ' mat2str(size(v)) ' ~= ' mat2str([n, 1])]);
            flag = 5;
            return
        elseif flag_direcao == 2
            flag = 5;
            return
        end
        if gerar_dados
            dados.vs = [dados.vs, v];
        end

        Jv = J * v;
        sigmaJv = opcoes.sigma * Jv;

        if gerar_dados
            dados.sigmaJvs = [dados.sigmaJvs, sigmaJv];
        end

        if opcoes.criterio_parada == 1
            criterio_parada = any(Jv >= -opcoes.epsilon);
        elseif opcoes.criterio_parada == 2
            criterio_parada = any(Jv + 0.5 * v' * v >= -opcoes.epsilon);
        end

        if criterio_parada
            if opcoes.mostra_cada_iteracao
                mostra_iteracao(xsol, fsol, k, [], v, sigmaJv, msg_v, toc(tic_iter_id), x0);
            end
            if gerar_dados
                dados.tempo = toc(tic_id);
            end
            if ~opcoes.quiet
                mostra(['Solução encontrada após ' num2str(k) ' iter. e ' num2str(toc(tic_id), '%.1f') ' seg.']);
            end
            flag = 1;
            return
        end

        % backtracking
        t = 1;
        while any(f(xsol + t * v) > fsol + t * sigmaJv) && t > 0
            t = t / 2;
            if gerar_dados
                dados.evalf = dados.evalf + 1;
            end
        end
        if gerar_dados
            dados.evalf = dados.evalf + 1;
            dados.ts = [dados.ts, t];
        end
        if opcoes.mostra_cada_iteracao
            mostra_iteracao(xsol, fsol, k, t, v, sigmaJv, msg_v, toc(tic_iter_id), x0);
        end
        if t == 0
            if gerar_dados
                dados.tempo = toc(tic_id);
            end
            mostra(['Backtracking com t = 0.']);
            flag = 5;
            return
        end

        xsol = xsol + t * v;
        fsol = f(xsol);
        if gerar_dados
            dados.fs = [dados.fs, fsol];
            dados.xs = [dados.xs, xsol];
        end
    end
    k = k + 1;
    if gerar_dados
        dados.tempo = toc(tic_id);
    end
    if ~opcoes.quiet
        mostra(['Máximo de iterações ou tempo máximo atingido, após ' num2str(k) ' iter. e ' num2str(toc(tic_id), '%.1f') ' seg.']);
    end
    if tem_callback
        opcoes.callback(xsol, k);
    end
    if opcoes.mostra_cada_iteracao
        mostra_iteracao(xsol, fsol, k, [], [], [], [], [], x0);
    end
    flag = 2;
end

function [] = mostra_iteracao(x, fx, k, t, v, sigmaJv, msg_v, tempo, x0)
    msg = [cor(37, 1) 'It. ' num2str(k) cor()];
    if ~isempty(v)
        if isempty(t)
            msg = [msg ' ' cor(32) '|v|= ' num2str(norm(v), 4) cor()];
            msg = [msg ' ' cor(33, 1) 'sigJv = ' vetor2str(sigmaJv, 4) cor()];
        else
            if t >= 1
                msg = [msg ' ' cor(31) 't= ' num2str(t) cor()];
            else
                msg = [msg ' ' cor(31) 't= 2^' num2str(log2(t)) cor()];
            end
            msg = [msg ' ' cor(32) '|tv|= ' num2str(norm(t * v), 4) cor()];
            msg = [msg ' ' cor(33, 1) 'tsigJv = ' vetor2str(t * sigmaJv, 4) cor()];
        end
        if ~isempty(msg_v)
            msg = [msg ' ' cor(33) msg_v cor()];
        end
    end
    if ~isempty(tempo)
        msg = [msg ' ' cor(34, 1) num2str(tempo, '%.1f') 's' cor()];
    end
    msg = [msg ' ' cor(35, 1) num2str(min(x), 3) ' <x< ' num2str(max(x), 3) cor()];
    msg = [msg ' ' cor(36) '|x - x0|= ' num2str(norm(x - x0), 3) cor()];
    msg = [msg ' ' cor(32, 1) 'f(x)= ' vetor2str(fx, 4) cor()];
    mostra(msg);
end
