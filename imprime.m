%
% imprime pontos em um arquivo para gerar gráfico.
%
% Entrada:
% xs       : pontos a serem impressos
% saida    : arquivo de saída
% undef    : (opcional, padrão []) pula linhas quando o ponto igual a undef
%

function [] = imprime(xs, saida, undef)
    % Constantes
    distancia_minima = 10^(-1.4);  % distância mínima entre pontos para serem plotados

    % Parâmetros
    usa_undef = nargin >= 3 && ~isempty(undef);

    saida_completa = ['gerado/' saida];
    if exist(saida_completa, 'file')
        mostra([cor(37, 2) 'Pulando ' saida_completa cor()]);
        return
    end
    if isempty(xs)
        mostra([cor(37, 2) 'Pulando arquivo vazio ' saida_completa cor()]);
        return
    end
    mostra([cor(37, 1) 'Gerando ' saida_completa cor()]);
    file = fopen(saida_completa, 'w');
    ponto_anterior = [];
    for k = 1:size(xs, 2)
        ponto_atual = xs(:, k);
        if ~usa_undef || any(ponto_atual ~= undef)
            linha = '';
            if isempty(ponto_anterior) || norm(ponto_anterior - ponto_atual) >= distancia_minima || k == size(xs, 2)
                for x = ponto_atual'
                    if ~isempty(linha)
                        linha = [linha ' '];
                    end
                    linha = [linha num2str(x)];
                end
                linha = [linha '\n'];
            end
        else
            linha = '\n';
        end
        if ~isempty(linha)
            fprintf(file, linha);
            ponto_anterior = ponto_atual;
        end
    end
    fclose(file);
end

