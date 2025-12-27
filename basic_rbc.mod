// ============================================
// RBC BÁSICO – MODELO MINIMAL
// ============================================

// Variables endógenas
var c k y a;

// Shock exógeno
varexo e;

// Parámetros
parameters alpha beta delta rho;

// Calibración simple
alpha = 0.36;
beta  = 0.99;
delta = 0.025;
rho   = 0.95;

// ============================================
// MODELO
// ============================================
model;
  // Producción
  y = a * k(-1)^alpha;

  // Restricción de recursos
  y = c + k - (1-delta)*k(-1);

  // Euler
  1/c = beta * (1/c(+1)) * (alpha*y(+1)/k + 1 - delta);

  // Proceso tecnológico
  a = rho*a(-1) + e;
end;

// ============================================
// VALORES INICIALES
// ============================================
initval;
  a = 1;
  k = 10;
  y = k^alpha;
  c = y - delta*k;
end;

// ============================================
// SHOCK
// ============================================
shocks;
  var e = 0.01^2;
end;

// ============================================
// RESOLUCIÓN
// ============================================
steady;
check;
stoch_simul(order=1, irf=20) y c k;
