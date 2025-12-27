// ============================================
// CAMPBELL (1994) - SECCIÓN 2
// Oferta laboral fija (L = 1)
// Versión pensada para replicar IRFs analíticas
// ============================================

// Variables endógenas (en logaritmos)
var y c k i a;

// Shock exógeno de tecnología (en logaritmos)
varexo e_a;

// Parámetros
parameters alpha beta delta rho_a sigma;

// Calibración estándar de Campbell (1994)
alpha = 0.667;   // participación del capital
beta  = 0.99;    // factor de descuento trimestral
delta = 0.025;   // depreciación trimestral
rho_a = 0.979;   // persistencia del shock tecnológico

// Elasticidad intertemporal de sustitución:
// sigma es el inverso de la EIS. Para log utilidad, sigma = 1.
sigma = 1;

// ============================================
// MODELO NO LINEAL
// ============================================
model;
  // (1) Producción: Y_t = A_t * K_{t-1}^alpha * L^(1-alpha), L=1
  exp(y) = exp(a) * exp(k(-1))^alpha;

  // (2) Restricción de recursos
  exp(y) = exp(c) + exp(i);

  // (3) Acumulación de capital
  exp(k) = (1-delta) * exp(k(-1)) + exp(i);

  // (4) Ecuación de Euler con CRRA, u(C) = C^(1-sigma)/(1-sigma)
  // u_c(C) = C^(-sigma)
  // C_t^(-sigma) = beta E_t[ C_{t+1}^(-sigma) * (r_{t+1} + 1 - delta) ]
  // r_{t+1} = alpha * Y_{t+1} / K_t
  exp(c)^(-sigma) = beta * exp(c(+1))^(-sigma) *
                    (alpha * exp(y(+1)) / exp(k) + 1 - delta);

  // (5) Proceso AR(1) para la tecnología
  a = rho_a * a(-1) + e_a;
end;

// ============================================
// ESTADO ESTACIONARIO ANALÍTICO
// ============================================
steady_state_model;
  // Tecnología en SS: a = 0  =>  A = 1
  a = 0;

  // 1. Capital de estado estacionario
  // Euler en estado estacionario:
  // 1 = beta * ( alpha * K_ss^(alpha-1) + 1 - delta )
  K_ss = ( alpha / ( 1/beta - (1-delta) ) )^(1/(1-alpha));
  k = log(K_ss);

  // 2. Producción
  Y_ss = K_ss^alpha;   // A=1, L=1
  y = log(Y_ss);

  // 3. Inversión: I_ss = delta * K_ss
  I_ss = delta * K_ss;
  i = log(I_ss);

  // 4. Consumo: C_ss = Y_ss - I_ss
  C_ss = Y_ss - I_ss;
  c = log(C_ss);
end;

// ============================================
// VALORES INICIALES
// ============================================
initval;
  a = 0;
  k = log( ( alpha / ( 1/beta - (1-delta) ) )^(1/(1-alpha)) );
  y = alpha * k;
  i = log(delta) + k;
  c = log(exp(y) - exp(i));
end;

resid;
steady;
check;

// ============================================
// SHOCK TECNOLÓGICO
// ============================================
shocks;
  var e_a = 0.0072^2; // misma varianza que en Campbell
end;

// ============================================
// IRFs Y SIMULACIÓN ESTOCÁSTICA
// ============================================
stoch_simul(order = 1, irf = 40) y c i k a;

// ============================================
// GRÁFICAS DE IRFs
// ============================================
// Gráficas automáticas de Dynare ya generadas arriba
// Gráficas personalizadas adicionales:
figure('Name', 'IRFs - Shock Tecnológico');
subplot(2,2,1); plot(0:40, [0; oo_.irfs.y_e_a(:)], 'b-', 'LineWidth', 2); title('Producto (y)'); grid on;
subplot(2,2,2); plot(0:40, [0; oo_.irfs.c_e_a(:)], 'r-', 'LineWidth', 2); title('Consumo (c)'); grid on;
subplot(2,2,3); plot(0:40, [0; oo_.irfs.i_e_a(:)], 'g-', 'LineWidth', 2); title('Inversión (i)'); grid on;
subplot(2,2,4); plot(0:40, [0; oo_.irfs.k_e_a(:)], 'm-', 'LineWidth', 2); title('Capital (k)'); grid on;

% Guardar figura en plots/
if ~exist('plots','dir')
  mkdir('plots');
end
print('-dpng','plots/campbel_sec2_irfs_tech.png');


