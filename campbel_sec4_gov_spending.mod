// ============================================
// CAMPBELL (1994) - SECCIÓN 4 - SOLUCIÓN ANALÍTICA
// ============================================

var y c k i n g;
varexo e_g;

parameters alpha beta delta sigma phi_n phi G_ss theta;

// PARÁMETROS
alpha = 0.667;
beta = 0.99;
delta = 0.025;
sigma = 1;          // Log utility
phi_n = 1.0;        // Elasticidad Frisch (ϕ)
phi = 0.95;         // Persistencia gasto (φ)
G_ss = 0.2;         // Gasto en SS (nivel)
// theta será calibrado en steady_state_model

model;
    // 1. PRODUCCIÓN
    exp(y) = exp(k(-1))^alpha * exp(n)^(1-alpha);

    // 2. RESTRICCIÓN RECURSOS
    // g es log(G/G_ss), entonces G_t = G_ss * exp(g)
    exp(y) = exp(c) + exp(i) + G_ss * exp(g);

    // 3. ACUMULACIÓN CAPITAL
    exp(k) = (1-delta) * exp(k(-1)) + exp(i);

    // 4. ECUACIÓN EULER
    exp(c)^(-sigma) = beta * exp(c(+1))^(-sigma) * 
                      (alpha * exp(y(+1))/exp(k) + 1 - delta);

    // 5. FOC TRABAJO CON theta
    theta * exp(n)^(1/phi_n) = (1/exp(c)) * (1-alpha) * exp(y) / exp(n);

    // 6. PROCESO GASTO
    g = phi * g(-1) + e_g;
end;

// ============================================
// ESTADO ESTACIONARIO ANALÍTICO
// ============================================
steady_state_model;
    // 1. RATIO Y/K DE EULER
    Y_over_K = (1/beta - (1-delta)) / alpha;
    
    // 2. NORMALIZACIÓN: N_ss = 0.33 (Campbell)
    N_ss = 0.33;
    n = log(N_ss);
    
    // 3. CAPITAL CONSISTENTE
    // De Y/K = ratio y Y = K^α * N^(1-α)
    K_ss = N_ss * (Y_over_K)^(1/(alpha-1));
    k = log(K_ss);
    
    // 4. PRODUCCIÓN
    Y_ss = K_ss^alpha * N_ss^(1-alpha);
    y = log(Y_ss);
    
    // 5. INVERSIÓN
    I_ss = delta * K_ss;
    i = log(I_ss);
    
    // 6. GASTO PÚBLICO (en desviaciones log)
    g = 0;  // log(G/G_ss) = 0 en SS
    
    // 7. CONSUMO DE RESTRICCIÓN
    C_ss = Y_ss - I_ss - G_ss;
    c = log(C_ss);
    
    // 8. CALIBRAR θ PARA QUE FOC TRABAJO SE CUMPLA
    // FOC: θ * N^(1/φ_n) = (1-α)*(Y/N) / C
    theta = ((1-alpha) * (Y_ss/N_ss) / C_ss) * (N_ss^(-1/phi_n));
end;

// ============================================
// VERIFICACIÓN
// ============================================
resid;
steady;
check;

// ============================================
// SHOCK
// ============================================
shocks;
    var e_g = 0.01^2;  // Shock 1% en gasto
end;

// ============================================
// SIMULACIÓN
// ============================================
stoch_simul(order = 1, irf = 40, nograph) y c i n g k;

// ============================================
// RESULTADOS
// ============================================
eta_yg = oo_.irfs.y_e_g(1);
eta_cg = oo_.irfs.c_e_g(1);
eta_ig = oo_.irfs.i_e_g(1);
eta_ng = oo_.irfs.n_e_g(1);
eta_kg = oo_.irfs.k_e_g(1);

printf('\n=============================================\n');
printf('RESULTADOS TABLA 6 - φ_n=%.1f, φ=%.2f\n', phi_n, phi);
printf('=============================================\n');
printf('Elasticidades impacto (t=0):\n');
printf('  η_{yx} (producto-gasto)   = %7.8f\n', eta_yg);
printf('  η_{cx} (consumo-gasto)    = %7.8f\n', eta_cg);
printf('  η_{ix} (inversión-gasto)  = %7.8f\n', eta_ig);
printf('  η_{nx} (trabajo-gasto)    = %7.8f\n', eta_ng);
printf('  η_{kx} (capital-gasto)    = %7.8f\n', eta_kg);
printf('=============================================\n');


// ============================================
// GRÁFICAS DE IRFs
// ============================================

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