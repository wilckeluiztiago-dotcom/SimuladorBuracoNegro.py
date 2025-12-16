// ============================================================
// SimuladorBuracoNegro - Constantes Físicas
// Autor: Luiz Tiago Wilcke
// ============================================================

#ifndef CONSTANTES_HPP
#define CONSTANTES_HPP

#include <cmath>

namespace BuracoNegro {

// ============================================================
// CONSTANTES FUNDAMENTAIS
// ============================================================

// Velocidade da luz no vácuo (m/s)
constexpr double VELOCIDADE_LUZ = 299792458.0;
constexpr double C = VELOCIDADE_LUZ;
constexpr double C2 = C * C;

// Constante gravitacional (m³/(kg·s²))
constexpr double CONSTANTE_GRAVITACIONAL = 6.67430e-11;
constexpr double G = CONSTANTE_GRAVITACIONAL;

// Constante de Planck (J·s)
constexpr double CONSTANTE_PLANCK = 6.62607015e-34;
constexpr double H_PLANCK = CONSTANTE_PLANCK;
constexpr double H_BARRA = H_PLANCK / (2.0 * M_PI);

// Constante de Boltzmann (J/K)
constexpr double CONSTANTE_BOLTZMANN = 1.380649e-23;
constexpr double K_BOLTZMANN = CONSTANTE_BOLTZMANN;

// Constante de Stefan-Boltzmann (W/(m²·K⁴))
constexpr double STEFAN_BOLTZMANN = 5.670374419e-8;

// Massa do Sol (kg)
constexpr double MASSA_SOL = 1.98892e30;

// Raio do Sol (m)
constexpr double RAIO_SOL = 6.96e8;

// Unidade astronômica (m)
constexpr double UNIDADE_ASTRONOMICA = 1.495978707e11;

// Parsec (m)
constexpr double PARSEC = 3.0856775814913673e16;

// Ano-luz (m)
constexpr double ANO_LUZ = 9.4607304725808e15;

// ============================================================
// CONSTANTES DERIVADAS PARA BURACOS NEGROS
// ============================================================

// Raio de Schwarzschild para 1 massa solar (m)
// rs = 2GM/c²
constexpr double RAIO_SCHWARZSCHILD_SOL = 2.0 * G * MASSA_SOL / C2;

// Fator para cálculo do raio de Schwarzschild
// rs = FATOR_RS * M (onde M em massas solares)
constexpr double FATOR_RAIO_SCHWARZSCHILD = 2.0 * G * MASSA_SOL / C2;

// Temperatura de Hawking para 1 massa solar (K)
// T = ℏc³ / (8πGMk)
constexpr double TEMP_HAWKING_SOL = H_BARRA * C * C * C / 
    (8.0 * M_PI * G * MASSA_SOL * K_BOLTZMANN);

// Luminosidade de Hawking para 1 massa solar (W)
// L = ℏc⁶ / (15360πG²M²)
constexpr double LUMINOSIDADE_HAWKING_SOL = H_BARRA * std::pow(C, 6) / 
    (15360.0 * M_PI * G * G * MASSA_SOL * MASSA_SOL);

// ============================================================
// FUNÇÕES UTILITÁRIAS INLINE
// ============================================================

// Raio de Schwarzschild para massa M (em kg)
inline double raio_schwarzschild(double massa_kg) {
    return 2.0 * G * massa_kg / C2;
}

// Raio de Schwarzschild para massa M (em massas solares)
inline double raio_schwarzschild_solar(double massa_solar) {
    return FATOR_RAIO_SCHWARZSCHILD * massa_solar;
}

// Temperatura de Hawking (K)
inline double temperatura_hawking(double massa_kg) {
    return H_BARRA * C * C * C / (8.0 * M_PI * G * massa_kg * K_BOLTZMANN);
}

// Luminosidade de Hawking (W)
inline double luminosidade_hawking(double massa_kg) {
    return H_BARRA * std::pow(C, 6) / 
           (15360.0 * M_PI * G * G * massa_kg * massa_kg);
}

// Tempo de evaporação (s)
inline double tempo_evaporacao(double massa_kg) {
    return 5120.0 * M_PI * G * G * std::pow(massa_kg, 3) / 
           (H_BARRA * std::pow(C, 4));
}

// Entropia de Bekenstein-Hawking
inline double entropia_bekenstein_hawking(double massa_kg) {
    double rs = raio_schwarzschild(massa_kg);
    double area = 4.0 * M_PI * rs * rs;
    return K_BOLTZMANN * C * C * C * area / (4.0 * G * H_BARRA);
}

// Raio do horizonte de eventos para Kerr (equatorial)
inline double raio_kerr(double massa_kg, double parametro_spin) {
    double M = G * massa_kg / C2;
    double a = parametro_spin * M;
    return M + std::sqrt(M * M - a * a);
}

// Ergosfera (raio externo no equador)
inline double raio_ergosfera(double massa_kg, double parametro_spin) {
    double M = G * massa_kg / C2;
    double a = parametro_spin * M;
    return M + std::sqrt(M * M - a * a);  // Simplificado para equador
}

// ISCO (Innermost Stable Circular Orbit) para Schwarzschild
inline double raio_isco_schwarzschild(double massa_kg) {
    return 3.0 * raio_schwarzschild(massa_kg);
}

// Raio da esfera de fótons (Schwarzschild)
inline double raio_esfera_fotons(double massa_kg) {
    return 1.5 * raio_schwarzschild(massa_kg);
}

// Gravidade superficial (Schwarzschild)
inline double gravidade_superficial(double massa_kg) {
    double rs = raio_schwarzschild(massa_kg);
    return C2 / (2.0 * rs);
}

} // namespace BuracoNegro

#endif // CONSTANTES_HPP
