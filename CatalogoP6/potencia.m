function archivo_potencia
  clc;
  #pkg load symbolic;
  warning('off', 'all');
  resultado = potencia([1.5 0.5; 0.5 1.5], [0; 1]);
endfunction

function resultado = potencia(A, x0)
    % Metodo de la potencia para obtener el autovalor dominante de una matriz A
    % :param A: Matriz A
    % :param x: x0 inicial
    % :param tol: Tolerancia al fallo que debe tener el resultado
    % :return: Autovalor dominante
    norma_inf_ant = 0;  % Resultado de la norma de la iteracion anterior
    itermax = 100;
    xk = x0;
    tol = 10**-10;
    errors = [];
    ks = [];
    for k=1:itermax
        yk = A * xk;  % Calculo de y
        % Calculo de la norma infinita
        ck = norm(yk, inf);
        % Calculo de xk
        xk_n = yk / ck;
        % Calculo del error
        error = norm(xk_n - xk);
        % Reasignacion de xk
        xk = xk_n
        errors = [errors error];
        ks = [ks k];
        % Tolerancia del error
        if error<tol 
            break;
        end
    end
    resultado = xk;
    plot(errors,ks, 'b')
    title("Grafica de puntos de comportamiento")
    ylabel("Iteraciones (k)")
    xlabel("Valores de los errores") 
end

