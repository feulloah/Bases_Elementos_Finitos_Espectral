function PHI = Base_Modal(P, alpha, beta)
%{
    Generación de Base de polinomios Modal, en el tramo [-1,1].

    Args:
        int  P : Valor del maximo grado de interpolacion
        float  alpha : Valor alpha de los polinomios de jacobi a utilizar.
                       Para base modal = 1. 
        float  beta : Valor beta de los polinomios de jacobi a utilizar.
                      Para base modal = 1.
    Returns:
        Array  PHI : Matriz de coeficientes de los polinomios Modales, 
                     donde cada fila representa un polinomio.
    %}

if nargin == 1
    alpha = 1;
    beta =  1;
end

PHI = zeros(P+1, P+1);
pjacs = Polinomios_Jacobi(P-1, alpha, beta);
for p = 0:P
    if p == 0
        PHI(p+1, end-1:end) = -.5*poly(1);
    elseif p == P
        PHI(p+1, end-1:end) = .5*poly(-1);
    else
        PHI(p+1, end-(p+1):end) = conv(conv(-.5*poly(1),.5*poly(-1)),pjacs(p, end-p+1:end));
    end
end