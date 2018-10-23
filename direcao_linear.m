%
% Encontra direção linear, resolvendo o problema
% argmin_(v,tau) tau
% s.a            J * v - tau <= 0
%                A * (x + v)  = a
%                B * (x + v) <= b
%             max(x - L, lx) <= x + v <= min(x + L, ux)
% onde L é a menor norma dos gradientes não nulos em J dividido pela raiz de n
%
% Entrada:
% implementa a função vf, parâmetro do metodo_descida
%
% Saída:
% implementa a função vf, parâmetro do metodo_descida
%

function [v, msg_v, flag] = direcao_linear(x, k, J, Jf, lx, ux, A, a, B, b, opcoes)
    [m, n] = size(J);

    L = norma_minima(J, opcoes.epsilon);
    if L <= opcoes.epsilon
        v = zeros(n, 1);
        msg_v = ['normin= ' num2str(L)];
        flag = 1;
        return
    end

    % Função a ser minimizada : c' * vtau
    % Restrição               : R * vtau <= r
    %                         : S * vtau  = s
    %                         : lb <= vtau <= ub
    c = sparse([zeros(n, 1); 1]);
    R = sparse([J, -ones(m, 1)]); % Jv - tau
    r = sparse(zeros(m, 1));      % 0
    if ~isempty(B)
        R = sparse([R; B, zeros(length(b), 1)]);  % Bv
        r = sparse([r; -B * x + b]);              % -Bx + b
    end
    S = [];
    s = [];
    if ~isempty(A)
        S = sparse([S; A, zeros(length(a), 1)]); % Av
        s = sparse([s; -A * x + a]);             % -Ax + a
    end
    if isempty(lx)
        lb = [-L * ones(n, 1); -999999999];
    else
        lb = [max(-L, lx - x); -999999999];
    end
    if isempty(ux)
        ub = [L * ones(n, 1); +999999999];
    else
        ub = [min(L, ux - x); +999999999];
    end
    RS = [R; S];
    rs = [r; s];
    ctype   = char(['U' * ones(size(R, 1), 1); 'S' * ones(size(S, 1), 1)]); % restrição do tipo <= e =
    vartype = char('C' * ones(n + 1, 1)); % variáveis contínuas
    sense   = 1; % indica que é problema de minimização

    % glpk é o método do Octave para resolver minimização linear com restrições lineares
    [vtau, fmin, errnum, extra] = glpk(c, RS, rs, lb, ub, ctype, vartype, sense);

    v     = vtau(1:n);
    tau   = vtau(n + 1);
    msg_v = ['tau/L= ' num2str(tau, 3) '/' num2str(L)];

    if errnum == 8
        mostra(['glpk atingiu o máximo de iterações permitidas']);
        flag = 2;
    elseif errnum ~= 0 || extra.status ~= 5
        mostra(['glpk não encontrou a solução. errnum = ' num2str(errnum) ' status = ' num2str(extra.status)]);
        flag = 2;
    else
        flag = 1;
    end
end
