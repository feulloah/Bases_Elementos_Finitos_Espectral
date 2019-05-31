function PHI = Base_NodalLagrange(P)
%{
    Generaci�n de Base de Polinomios de interpolaci�n de lagrange de grado 
    P con nodos equiespaciados en el tramo [-1, 1]

    Args:
        int  P : Grado m�ximo de la interpolaci�n. 
    Returns:
        Array  L : Matriz de coeficientes de los polinomios de Lagrange,
        donde cada fila representa un polinomio.
%}
PHI = zeros(P+1, P+1);
chis = linspace(-1,1,P+1);
for i = 1:P+1
    V = 1;
    for j = 1:P+1
        if j~=i
            V=conv(V,poly(chis(j)))/(chis(i)-chis(j));
        end
    end
    PHI(i, :) = V;
end
