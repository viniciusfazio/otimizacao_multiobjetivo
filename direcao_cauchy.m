%
% Encontra direção de Cauchy, resolvendo o problema
% argmin_(v,tau) 0.5 * norm(v) + tau
% s.a            J * v - tau <= 0
%                A * (x + v)  = a
%                B * (x + v) <= b
%                lx <= x + v <= ux
%
% Entrada:
% implementa a função vf, parâmetro do metodo_descida
%
% Saída:
% implementa a função vf, parâmetro do metodo_descida
%

function [v, msg_v, flag] = direcao_cauchy(x, k, J, Jf, lx, ux, A, a, B, b, opcoes)
    [m, n] = size(J);

    % Se todos gradientes forem factíveis e m <= 2, resolve o problema de forma mais direta
    %mostra('verificando forma direta');
    JsatisfazA = isempty(a) || all(all(abs(A * -J' + repmat(A * x - a, 1, m)) <= opcoes.epsilon));
    JsatisfazB = isempty(b) || all(all(B * -J' + repmat(B * x - b, 1, m) <= 0));
    Jsatisfazlx = isempty(lx) || all(all(-J' + repmat(x - lx, 1, m) >= 0));
    Jsatisfazux = isempty(ux) || all(all(-J' + repmat(x - ux, 1, m) <= 0));
    if m <= 2 && JsatisfazA && JsatisfazB && Jsatisfazlx && Jsatisfazux
        flag = 1;
        msg_v = 'direto';
        if m == 1
            v = -J';
            return
        end
        j1 = J(1, :)';
        j2 = J(2, :)';
        d = j1 - j2;
        if norm(d) < opcoes.epsilon
            v = -j1;
        else
            c = j1' * d / (d' * d);
            if c <= 0
                v = -j1;
            elseif c >= 1
                v = -j2;
            else
                v = -j1 + c * d ;
            end
        end
        return
    end

    % Função a ser minimizada : 0.5vtau' * H * vtau + vtau' * q
    % Restrição               : R * vtau <= r
    %                         : S * vtau  = s
    %                         : lb <= vtau <= ub
    % Ponto inicial           : vtau0 = 0
    H = sparse(diag([ones(n, 1) ; 0]));
    q = sparse([zeros(n, 1); 1]);
    R = sparse([J, -ones(m, 1)]); % Jv - tau
    r = sparse(zeros(m, 1));      % 0
    if ~isempty(B)
        R = sparse([R; B, zeros(length(b), 1)]); % Bv
        r = sparse([r; -B * x + b]);             % -Bx + b
    end
    S = [];
    s = [];
    if ~isempty(A)
        S = sparse([S; A, zeros(length(a), 1)]); % Av
        s = sparse([s; -A * x + a]);             % -Ax + a
    end
    lb = [];
    if ~isempty(lx)
        lb = [lb; lx - x; -999999999];
    end
    ub = [];
    if ~isempty(ux)
        ub = [ub; ux - x; +999999999];
    end
    vtau0 = sparse(zeros(n + 1, 1));

    % qp é o método do Octave para resolver minimização quadrática com restrições lineares
    [vtau, obj, info] = qp(vtau0, H, q, S, s, lb, ub, [], R, r);

    v     = vtau(1:n);
    tau   = vtau(n + 1);
    msg_v = ['tau= ' num2str(tau, 3)];

    flag = 1;
    if info.info == 3
        mostra(['qp atingiu o máximo de iterações permitidas']);
        flag = 2;
    elseif info.info ~= 0
        mostra(['qp não encontrou a solução. info = ' num2str(info.info)]);
        flag = 2;
    end
    if flag ~= 1
        satisfazA = isempty(a) || all(abs(A * (x + v) - a) <= opcoes.epsilon);
        satisfazB = isempty(b) || all(B * (x + v) <= b);
        satisfazlx = isempty(lx) || all(x + v >= lx);
        satisfazux = isempty(ux) || all(x + v <= ux);
        if opcoes.criterio_parada == 1
            if any(J * v >= -opcoes.epsilon) || ~satisfazA || ~satisfazB || ~satisfazlx || ~satisfazux || any(isna(v))
                return;
            end
        elseif opcoes.criterio_parada == 2
            if any(J * v + 0.5 * v' * v >= -opcoes.epsilon) || ~satisfazA || ~satisfazB || ~satisfazlx || ~satisfazux || any(isna(v))
                return;
            end
        end
        mostra('erro ignorado.');
        flag = 1;
    end
end
