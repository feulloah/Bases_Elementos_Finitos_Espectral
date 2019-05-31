function PHI = Base_NodalEspectral(p)
    %{
    Generación de Base de polinomios Nodal Espectral, en el tramo [-1,1].

    Args:
        float  p : Grado del polinomio de interpolación.
    Returns:
        Array  PHI : Matriz de coeficientes de los polinomios de Nodal 
                    Espectral, donde cada fila representa un polinomio.
    %}
PHI = zeros(p+1, p+1);
Lp = Polinomios_Jacobi(p, 0, 0);
Lp = Lp(end, :);
dL = zeros(size(Lp));
dL(2:end)=polyder(Lp);
g_p = Polinomios_Jacobi(p, 1, 1);
g_p = g_p(end-1, :);
g_p = conv(conv(g_p, poly(1)),poly(-1));
chi_p = sort(roots(g_p));
for i = 1:p+1
    hp_0 = conv(conv(poly(-1),poly(1)), dL(2:end));
    [PHI(i, :), r] = deconv(hp_0, ((p*(p+1)*polyval(Lp, chi_p(i)))*poly(chi_p(i))));
    r(1) = [];
    PHI(i, :) = PHI(i, :) + r;
end