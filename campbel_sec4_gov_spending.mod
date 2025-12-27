// ============================================
// CAMPBELL (1994) - SECCIÓN 4
// Choques de gasto público con oferta laboral variable
// Dos esquemas de financiación:
//   1) Lump-sum (tau = 0)  -> no distorsiona FOCs, solo RC
//   2) Impuesto proporcional al producto (tau > 0)
// ============================================

var y c k i n a g;     // todas en logaritmos
varexo e_a e_g;        // shocks de tecnología y gasto público

parameters alpha beta delta rho_a rho_g sigma phi_n tau g_y;

// Tecnología y preferencias
alpha = 0.667;
beta  = 0.99;
delta = 0.025;
rho_a = 0.979;
rho_g = 0.9;

sigma = 1;        // inverso de EIS (log utilidad en c)
phi_n = 1.0;      // peso del ocio (separable)

// Gasto público en estado estacionario (como fracción de Y)
g_y = 0.2;

// Tipo impositivo proporcional sobre el producto:
//   tau = 0   -> financiación lump-sum (sin distorsión)
//   tau > 0   -> financiación distorsionante
tau = 0.0;

// ============================================
// MODELO NO LINEAL
// ============================================
model;
  // (1) Producción
  exp(y) = exp(a) * exp(k(-1))^alpha * exp(n)^(1-alpha);

  // (2) Restricción de recursos agregada:
  // Y_t = C_t + I_t + G_t
  exp(y) = exp(c) + exp(i) + exp(g);

  // (3) Acumulación de capital
  exp(k) = (1-delta) * exp(k(-1)) + exp(i);

  // (4) Euler intertemporal con posible impuesto distorsionante
  // u_c(C_t) = C_t^(-sigma)
  // u_c(C_t) = beta E_t[ u_c(C_{t+1}) * ( (1-tau)*r_{t+1} + 1 - delta ) ]
  // r_{t+1} = alpha * Y_{t+1} / K_t
  exp(c)^(-sigma) = beta * exp(c(+1))^(-sigma) *
                    ( (1-tau) * alpha * exp(y(+1)) / exp(k) + 1 - delta );

  // (5) FOC intratemporal trabajo-ocio (modelo separable)
  // u(c,n) = log c + phi_n log(1-n)
  // phi_n / (1-n) = (1/c) * (1-tau) * w_t
  // w_t = (1-alpha) * Y_t / N_t
  phi_n / (1 - exp(n)) = (1/exp(c)) * (1-tau) * (1-alpha) * exp(y) / exp(n);

  // (6) Proceso de tecnología
  a = rho_a * a(-1) + e_a;

  // (7) Proceso de gasto público (en log de desviaciones respecto al ss)
  // g_hat_t = rho_g * g_hat_{t-1} + e_g
  g = rho_g * g(-1) + e_g;
end;

// ============================================
// ESTADO ESTACIONARIO (Dynare lo resolverá numéricamente)
// ============================================
// Nota: El estado estacionario se resuelve numéricamente
// porque n, k y g están acoplados a través de las FOCs

initval;
  a = 0;
  g = 0;
  n = log(0.33);
  k = log( ( alpha / ( 1/beta - (1-delta) ) )^(1/(1-alpha)) );
  y = log( exp(k)^alpha * exp(n)^(1-alpha) );
  i = log(delta) + k;
  c = log(exp(y) - exp(i) - g_y * exp(y));
end;

resid;
steady;
check;

shocks;
  var e_a = 0.0072^2;
  var e_g = 0.01^2;
end;

// ============================================
// IRFs Y SIMULACIÓN ESTOCÁSTICA
// ============================================
stoch_simul(order = 1, irf = 40) y c i k n a g;

// ============================================
// GRÁFICAS DE IRFs
// ============================================
// Gráficas automáticas de Dynare ya generadas arriba
// IRFs a shock tecnológico
figure('Name', 'IRFs - Shock Tecnológico');
subplot(2,3,1); plot(0:40, [0; oo_.irfs.y_e_a(:)], 'b-', 'LineWidth', 2); title('Producto (y)'); grid on;
subplot(2,3,2); plot(0:40, [0; oo_.irfs.c_e_a(:)], 'r-', 'LineWidth', 2); title('Consumo (c)'); grid on;
subplot(2,3,3); plot(0:40, [0; oo_.irfs.i_e_a(:)], 'g-', 'LineWidth', 2); title('Inversión (i)'); grid on;
subplot(2,3,4); plot(0:40, [0; oo_.irfs.k_e_a(:)], 'm-', 'LineWidth', 2); title('Capital (k)'); grid on;
subplot(2,3,5); plot(0:40, [0; oo_.irfs.n_e_a(:)], 'c-', 'LineWidth', 2); title('Trabajo (n)'); grid on;
subplot(2,3,6); plot(0:40, [0; oo_.irfs.a_e_a(:)], 'k-', 'LineWidth', 2); title('Tecnología (a)'); grid on;

% Guardar figura de shock tecnológico
if ~exist('plots','dir')
  mkdir('plots');
end
print('-dpng','plots/campbel_sec4_irfs_tech.png');

// IRFs a shock de gasto público
figure('Name', 'IRFs - Shock de Gasto Público');
subplot(2,3,1); plot(0:40, [0; oo_.irfs.y_e_g(:)], 'b-', 'LineWidth', 2); title('Producto (y)'); grid on;
subplot(2,3,2); plot(0:40, [0; oo_.irfs.c_e_g(:)], 'r-', 'LineWidth', 2); title('Consumo (c)'); grid on;
subplot(2,3,3); plot(0:40, [0; oo_.irfs.i_e_g(:)], 'g-', 'LineWidth', 2); title('Inversión (i)'); grid on;
subplot(2,3,4); plot(0:40, [0; oo_.irfs.k_e_g(:)], 'm-', 'LineWidth', 2); title('Capital (k)'); grid on;
subplot(2,3,5); plot(0:40, [0; oo_.irfs.n_e_g(:)], 'c-', 'LineWidth', 2); title('Trabajo (n)'); grid on;
subplot(2,3,6); plot(0:40, [0; oo_.irfs.g_e_g(:)], 'k-', 'LineWidth', 2); title('Gasto Público (g)'); grid on;

% Guardar figura de shock fiscal
if ~exist('plots','dir')
  mkdir('plots');
end
print('-dpng','plots/campbel_sec4_irfs_fiscal.png');


