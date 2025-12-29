// ============================================
// CAMPBELL (1994) - SECCIÓN 4 
// ============================================

var y c k i n g;
varexo e_g;

parameters alpha beta delta sigma phi_n phi G_ss theta;

// PARÁMETROS EXACTOS DE CAMPBELL
alpha = 0.667;      // α = participación TRABAJO (p. 467)
beta = 0.99;        // Factor descuento trimestral
delta = 0.025;      // Depreciación trimestral
sigma = 1;          // Elasticidad intertemporal sustitución (log utility)
phi_n = 1.0;        // Elasticidad Frisch (ϕ)
phi = 0.95;         // Persistencia gasto (φ)
G_ss = 0.2;         // Gasto en SS (20% del PIB)

model;
    // 1. PRODUCCIÓN (Campbell eq. 1, p. 467)
    // Y_t = N_t^α * K_{t-1}^{1-α}
    exp(y) = exp(n)^alpha * exp(k(-1))^(1-alpha);

    // 2. RESTRICCIÓN RECURSOS (Campbell eq. 50, p. 480)
    // Y_t = C_t + I_t + G_t
    // g_t = log(G_t/G_ss) → G_t = G_ss * exp(g_t)
    exp(y) = exp(c) + exp(i) + G_ss * exp(g);

    // 3. ACUMULACIÓN CAPITAL (Campbell eq. 2, p. 467)
    exp(k) = (1-delta) * exp(k(-1)) + exp(i);

    // 4. ECUACIÓN EULER (Campbell eq. 36, p. 475 con τ=0)
    // MPK = (1-α) * Y_{t+1}/K_t
    exp(c)^(-sigma) = beta * exp(c(+1))^(-sigma) * 
                      ((1-alpha) * exp(y(+1))/exp(k) + 1 - delta);

    // 5. FOC TRABAJO (Campbell eq. 36, p. 475 con τ=0)
    // MPL = α * Y_t/N_t
    // U(C,N) = log(C) - θ * N^{1+1/ϕ}/(1+1/ϕ)
    theta * exp(n)^(1/phi_n) = (1/exp(c)) * alpha * exp(y) / exp(n);

    // 6. PROCESO GASTO (Campbell p. 481)
    // log(G_t) = φ * log(G_{t-1}) + (1-φ)*log(G_ss) + ε_t^g
    // En desviaciones: g_t = φ * g_{t-1} + e_g
    g = phi * g(-1) + e_g;
end;

// ============================================
// ESTADO ESTACIONARIO ANALÍTICO
// ============================================
steady_state_model;
    // 1. RATIO Y/K DE EULER (Campbell eq. 10, p. 468)
    // MPK = (1-α)Y/K = 1/β - (1-δ)
    Y_over_K = (1/beta - (1-delta)) / (1-alpha);
    
    // 2. NORMALIZACIÓN: N_ss = 0.33 (1/3 tiempo trabajando)
    N_ss = 0.33;
    n = log(N_ss);
    
    // 3. CAPITAL CONSISTENTE
    // De Y = N^α * K^{1-α} y Y/K = ratio
    // K^{α-1} * N^α = Y/K
    // K = N * (Y/K)^{1/(α-1)}
    K_ss = N_ss * (Y_over_K)^(1/(alpha-1));
    k = log(K_ss);
    
    // 4. PRODUCCIÓN
    Y_ss = N_ss^alpha * K_ss^(1-alpha);
    y = log(Y_ss);
    
    // 5. INVERSIÓN
    I_ss = delta * K_ss;
    i = log(I_ss);
    
    // 6. GASTO PÚBLICO (en desviaciones log)
    g = 0;  // log(G/G_ss) = 0 en SS
    
    // 7. CONSUMO (de restricción recursos)
    C_ss = Y_ss - I_ss - G_ss;
    c = log(C_ss);
    
    // 8. CALIBRAR θ PARA FOC TRABAJO (Campbell eq. 36)
    // θ * N^(1/ϕ) = α * (Y/N) / C
    theta = (alpha * (Y_ss/N_ss) / C_ss) * (N_ss^(-1/phi_n));
end;

// ============================================
// VERIFICACIÓN
// ============================================
resid;
steady;
check;

// ============================================
// SHOCK (1% en gasto)
// ============================================
shocks;
    var e_g = 0.01^2;  // Shock 1% en log(G)
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
printf('CAMPBELL (1994) - TABLA 6 - NOTACIÓN EXACTA\n');
printf('α=%.3f (trabajo), 1-α=%.3f (capital)\n', alpha, 1-alpha);
printf('ϕ=%.1f (Frisch), φ=%.2f (persistencia)\n', phi_n, phi);
printf('=============================================\n');
printf('ESTADO ESTACIONARIO:\n');
printf('  Y_ss = %.4f,  K_ss = %.4f,  N_ss = %.4f\n', exp(y), exp(k), exp(n));
printf('  C_ss = %.4f,  I_ss = %.4f,  G_ss = %.4f\n', exp(c), exp(i), G_ss);
printf('  C/Y = %.3f,  I/Y = %.3f,  G/Y = %.3f\n', exp(c)/exp(y), exp(i)/exp(y), G_ss/exp(y));
printf('  θ = %.4f,  Y/K = %.4f\n', theta, exp(y)/exp(k));
printf('\nELASTICIDADES IMPACTO (η_{yx}, etc.):\n');
printf('  Producto:    %7.8f\n', eta_yg);
printf('  Consumo:     %7.8f\n', eta_cg);
printf('  Inversión:   %7.8f\n', eta_ig);
printf('  Trabajo:     %7.8f\n', eta_ng);
printf('  Capital:     %7.8f\n', eta_kg);
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
