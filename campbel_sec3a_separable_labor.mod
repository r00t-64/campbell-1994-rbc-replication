// ============================================
// CAMPBELL (1994) - SECCIÓN 3a
// Oferta laboral variable, utilidad separable:
// u(c, n) = log c + phi_n * log(1-n)
// ============================================

var y c k i n a;      // todas en logaritmos salvo que se indique
varexo e_a;           // shock tecnológico

parameters alpha beta delta rho_a sigma phi_n;

// Calibración base (puede ajustarse a la del paper)
alpha = 0.667;
beta  = 0.99;
delta = 0.025;
rho_a = 0.979;

// Inverso de la elasticidad intertemporal de sustitución
sigma = 1;           // log utilidad en consumo

// Parámetro de utilidad del ocio (relacionado con elasticidad de la oferta laboral)
phi_n = 1.0;

// ============================================
// MODELO NO LINEAL
// ============================================
model;
  // (1) Producción: Y_t = A_t * K_{t-1}^alpha * N_t^(1-alpha)
  exp(y) = exp(a) * exp(k(-1))^alpha * exp(n)^(1-alpha);

  // (2) Restricción de recursos
  exp(y) = exp(c) + exp(i);

  // (3) Acumulación de capital
  exp(k) = (1-delta) * exp(k(-1)) + exp(i);

  // (4) Ecuación de Euler intertemporal (CRRA en c)
  // C_t^(-sigma) = beta E_t[ C_{t+1}^(-sigma) * (r_{t+1} + 1 - delta) ]
  // r_{t+1} = alpha * Y_{t+1} / K_t
  exp(c)^(-sigma) = beta * exp(c(+1))^(-sigma) *
                    (alpha * exp(y(+1)) / exp(k) + 1 - delta);

  // (5) FOC estática trabajo-ocio (modelo separable)
  // u(c,n) = log c + phi_n * log(1-n)
  // u_c = 1/c,  u_n = -phi_n / (1-n)
  // Condición: -u_n / u_c = w_t
  //  => phi_n / (1-n) = (1/c) * w_t
  // Salario real: w_t = (1-alpha) * Y_t / N_t
  phi_n / (1 - exp(n)) = (1/exp(c)) * (1-alpha) * exp(y) / exp(n);

  // (6) Proceso AR(1) para la tecnología
  a = rho_a * a(-1) + e_a;
end;

// ============================================
// ESTADO ESTACIONARIO (Dynare lo resolverá numéricamente)
// ============================================
// Nota: El estado estacionario se resuelve numéricamente
// porque n y k están acoplados a través de las FOCs

initval;
  a = 0;
  n = log(0.33);
  k = log( ( alpha / ( 1/beta - (1-delta) ) )^(1/(1-alpha)) );
  y = log( exp(k)^alpha * exp(n)^(1-alpha) );
  i = log(delta) + k;
  c = log(exp(y) - exp(i));
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
figure('Name', 'IRFs - Shock Tecnológico (Oferta Laboral Variable)');
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
print('-dpng','plots/campbel_sec3a_irfs_tech.png');


