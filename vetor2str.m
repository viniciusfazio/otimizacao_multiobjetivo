%
% Transforma um vetor m em texto, resumindo o conteúdo se a dimensão for grande
%
% Entrada:
% m        : vetor coluna
% precisao : número de dígitos significativos
%
% Saída:
% s : texto
%

function s = vetor2str(m, precisao)
    if size(m, 1) > 5
        s = ['mi.me.ma_' mat2str([min(m) mean(m) max(m)], precisao)];
    else
        s = mat2str(m', precisao);
    end
end
