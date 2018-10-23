%
% Retorna a string para mudar a cor do texto para [n1;n2m na saída Bash
%

function s = cor(n1, n2)
    if nargin < 1
        s = "\033[0m";
        return
    end
    if nargin < 2
        n2 = 0;
    end
    s = ["\033[" num2str(n2) ";" num2str(n1) "m"];
end
