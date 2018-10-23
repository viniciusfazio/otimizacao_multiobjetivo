%
% Imprime um texto e faz flush para mostrar imediatamente
%

function [] = mostra(msg)
    if nargin < 1
        msg = '';
    end
    disp(msg);
    fflush(stdout);
end
