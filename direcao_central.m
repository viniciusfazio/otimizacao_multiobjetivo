%
% Encontra direção de Cauchy, resolvendo o problema
% argmin_(v,tau) 0.5 * norm(v) + tau
% s.a            DJ * v - tau <= 0
%                 A * (x + v)  = a
%                 B * (x + v) <= b
%                 lx <= x + v <= ux
% onde D é uma matriz diagonal positiva que normaliza os gradientes de J pela menor norma entre os gradientes
%
% Entrada:
% implementa a função vf, parâmetro do metodo_descida
%
% Saída:
% implementa a função vf, parâmetro do metodo_descida
%

function [v, msg_v, flag] = direcao_central(x, k, J, Jf, lx, ux, A, a, B, b, opcoes)
    [L, DJ] = norma_minima(J, opcoes.epsilon);
    if L <= opcoes.epsilon
        v = x * 0;
        msg_v = ['normin= ' num2str(L)];
        flag = 1;
        return
    end

    [v, msg_v, flag] = direcao_cauchy(x, k, DJ, [], lx, ux, A, a, B, b, opcoes);
end
