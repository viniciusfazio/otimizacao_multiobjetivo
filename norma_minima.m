%
% Retorna a menor norma entre os gradientes, e opcionalmente o J normalizado por esta norma
%
% Entrada:
% J       : matriz jacobiana
% epsilon : valor da norma mínimo para que seja feita a normalização
%
% Saida:
% menor_norma  : valor da norma mínima
% Jnormalizado : matriz jacobiana normalizada pela norma mínima

function [menor_norma, Jnormalizado] = norma_minima(J, epsilon)
    [m, n] = size(J);

    menor_norma = -1;
    if nargout > 1
        Jnormalizado = zeros(m, n);
    end
    for k = 1:m
        norma = norm(J(k, :));
        if menor_norma > norma || menor_norma < 0
            menor_norma = norma;
        end
        if menor_norma <= epsilon
            return
        end
        if nargout > 1
            Jnormalizado(k, :) = J(k, :) / norma;
        end
    end
    if nargout > 1
        Jnormalizado = Jnormalizado * menor_norma;
    end
end
