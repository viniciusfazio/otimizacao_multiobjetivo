%
% Encontra direção de Pareto descendente, resolvendo o problema
% argmin_(w,tau) tau
% s.a            DJ * DJ' * w + tau >= 0
%                            1' * w  = 1
%                 A * (x - DJ' * w)  = a
%                 B * (x - DJ' * w) <= b
%                 lx <= x - DJ' * w <= ux
%                                 w >= 0
% onde D é uma matriz diagonal positiva que normaliza os gradientes de J pela menor norma entre os gradientes
% e 1' é um vetor linha de m componentes iguais a 1.
%
% Se tau > -beta * L, retorna direção linear, onde M é a norma do menor gradiente e beta é uma constante entre 0 e 1.
%
% Entrada:
% implementa a função vf, parâmetro do metodo_descida
% opcoes.beta : (padrão 0.1) parâmetro da "condição do ângulo" da direção de Pareto descendente
%
% Saída:
% implementa a função vf, parâmetro do metodo_descida


function [v, msg_v, flag] = direcao_pareto(x, k, J, Jf, lx, ux, A, a, B, b, opcoes)
    if ~isfield(opcoes, 'beta')
        opcoes.beta = 1e-1;
    end

    [m, n] = size(J);

    [L, DJ] = norma_minima(J, opcoes.epsilon);
    if L <= opcoes.epsilon
        flag = 1;
        v = zeros(n, 1);
        msg_v = ['normin= ' num2str(L)];
        return
    end

    % Função a ser minimizada : c' * wtau
    % Restrição               : R * wtau <= r
    %                         : S * wtau  = s
    %                         : lb <= wtau
    c = sparse([zeros(m, 1); 1]);
    R = sparse([-DJ * DJ', -ones(m, 1)]); % -DJDJ'w - tau
    r = sparse(zeros(m, 1));              % 0
    R = sparse([R; ones(1, m), 0]); % 1'w
    r = sparse([r; 1]);             % 1
    if ~isempty(B)
        R = sparse([R; -B * DJ', zeros(length(b), 1)]); % -BDJ'w
        r = sparse([r; -B * x + b]);                    % -Bx + b
    end
    if ~isempty(lx)
        R = sparse([R; DJ', zeros(n, 1)]); % DJ'w
        r = sparse([r; x - lx]);           % x - lx
    end
    if ~isempty(ux)
        R = sparse([R; -DJ', zeros(n, 1)]); % -DJ'w
        r = sparse([r; ux - x]);            % ux - x
    end
    S = [];
    s = [];
    if ~isempty(A)
        S = sparse([S; -A * DJ', zeros(length(a), 1)]); % -ADJ'w
        s = sparse([s; -A * x + a]);                    % -Ax + a
    end
    lb = [zeros(m, 1); -999999999];
    RS = [R; S];
    rs = [r; s];
    ctype   = char(['U' * ones(size(R, 1), 1); 'S' * ones(size(S, 1), 1)]); % restrição do tipo <= e =
    vartype = char('C' * ones(m + 1, 1)); % variáveis contínuas
    sense   = 1; % indica que é problema de minimização

    % glpk é o método do Octave para resolver minimização linear com restrições lineares
    [wtau, fmin, errnum, extra] = glpk(c, RS, rs, lb, [], ctype, vartype, sense);

    w   = wtau(1:m);
    tau = wtau(m + 1);
    if errnum ~= 0 || extra.status ~= 5 || tau > -opcoes.beta * L * L
        [v, msg_v, flag] = direcao_linear(x, k, J, Jf, lx, ux, A, a, B, b, opcoes);
        msg_v = ['(L)' msg_v];
    else
        msg_v = ['tau ' num2str(tau, 3)];
        v = -DJ' * w;
        flag = 1;
    end
end
