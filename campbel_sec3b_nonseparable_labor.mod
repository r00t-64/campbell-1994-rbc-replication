// ============================================
// CAMPBELL (1994) - SECCIÓN 3b
// Oferta laboral variable, utilidad no separable:
// u(c, n) = ( c^gamma * (1-n)^(1-gamma) )^(1-sigma) / (1-sigma)
// ============================================

var y c k i n a;
varexo e_a;

parameters alpha beta delta rho_a sigma gamma_u;

alpha = 0.667;
beta  = 0.99;
delta = 0.025;
rho_a = 0.979;

// Parámetros de preferencia
sigma   = 1;      // inverso de EIS
gamma_u = 0.5;    // peso del consumo en la utilidad compuesta

// ============================================
// MODELO NO LINEAL
// ============================================
model;
  // (1) Producción
  exp(y) = exp(a) * exp(k(-1))^alpha * exp(n)^(1-alpha);

  // (2) Restricción de recursos
  exp(y) = exp(c) + exp(i);

  // (3) Acumulación de capital
  exp(k) = (1-delta) * exp(k(-1)) + exp(i);

  // (4) Euler intertemporal (mismo que con utilidad separable en c)
  exp(c)^(-sigma) = beta * exp(c(+1))^(-sigma) *
                    (alpha * exp(y(+1)) / exp(k) + 1 - delta);

  // (5) FOC intratemporal con utilidad Cobb-Douglas no separable
  // u(c,n) = [c^gamma_u (1-n)^(1-gamma_u)]^(1-sigma)/(1-sigma)
  // Cociente marginal: -u_n / u_c = ((1-gamma_u)/gamma_u) * c / (1-n)
  // Condición: -u_n/u_c = w_t
  // w_t = (1-alpha) * Y_t / N_t
  ((1-gamma_u)/gamma_u) * exp(c) / (1 - exp(n))
    = (1-alpha) * exp(y) / exp(n);

  // (6) Proceso de tecnología
  a = rho_a * a(-1) + e_a;
end;

// ============================================
// ESTADO ESTACIONARIO ANALÍTICO
// Se obtiene una solución cerrada usando la condición de Euler en SS y la FOC intratemporal
steady_state_model;
  // Tecnología en SS
  a = 0;

  // kappa = Y/K en estado estacionario (de Euler)
  kappa = (1/beta - (1-delta)) / alpha;

  // Para utilidad no separable con peso gamma_u, la FOC intratemporal en SS da:
  // ((1-gamma_u)/gamma_u) * C/(1-N) = (1-alpha) * Y / N
  // con C = (kappa - delta) * K y Y = kappa * K, se obtiene:
  // N = B / (A + B) con
  // B = (1-alpha)*kappa
  // A = ((1-gamma_u)/gamma_u) * (kappa - delta)
  B = (1-alpha) * kappa;
  A = ((1 - gamma_u) / gamma_u) * (kappa - delta);
  N_ss = B / (A + B);
  n = log(N_ss);

  // Capital en SS
  K_ss = N_ss * kappa^(-1/(1-alpha));
  k = log(K_ss);

  // Producción, inversión y consumo en SS
  Y_ss = K_ss^alpha * N_ss^(1-alpha);
  y = log(Y_ss);
  i = log(delta) + k;
  c = log(Y_ss - exp(i));
end;

resid;
steady;
check;

shocks;
  var e_a = 0.0072^2;
end;

stoch_simul(order = 1, irf = 40) y c i k n a;

// ============================================
// GRÁFICAS DE IRFs
// ============================================
// Gráficas automáticas de Dynare ya generadas arriba
// Gráficas personalizadas adicionales:
figure('Name', 'IRFs - Shock Tecnológico (Utilidad No Separable)');
subplot(2,3,1); plot(0:40, [0; oo_.irfs.y_e_a(:)], 'b-', 'LineWidth', 2); title('Producto (y)'); grid on;
subplot(2,3,2); plot(0:40, [0; oo_.irfs.c_e_a(:)], 'r-', 'LineWidth', 2); title('Consumo (c)'); grid on;
subplot(2,3,3); plot(0:40, [0; oo_.irfs.i_e_a(:)], 'g-', 'LineWidth', 2); title('Inversión (i)'); grid on;
subplot(2,3,4); plot(0:40, [0; oo_.irfs.k_e_a(:)], 'm-', 'LineWidth', 2); title('Capital (k)'); grid on;
subplot(2,3,5); plot(0:40, [0; oo_.irfs.n_e_a(:)], 'c-', 'LineWidth', 2); title('Trabajo (n)'); grid on;
subplot(2,3,6); plot(0:40, [0; oo_.irfs.a_e_a(:)], 'k-', 'LineWidth', 2); title('Tecnología (a)'); grid on;

% Guardar figura
if ~exist('plots','dir')
  mkdir('plots');
end
print('-dpng','plots/campbel_sec3b_irfs_tech.png');


