// ============================================================
// SimuladorBuracoNegro - Métrica de Kerr
// Autor: Luiz Tiago Wilcke
// 
// Buraco negro rotativo (com momento angular)
// ============================================================

#ifndef KERR_HPP
#define KERR_HPP

#include "constantes.hpp"
#include <cmath>
#include <array>

namespace BuracoNegro {

class MetricaKerr {
private:
    double M_;      // Massa (unidades geométricas)
    double a_;      // Parâmetro de spin a = J/(Mc)
    double massa_kg_;
    double spin_;   // Parâmetro adimensional a/M

    // Funções auxiliares
    double Sigma(double r, double theta) const {
        return r * r + a_ * a_ * std::cos(theta) * std::cos(theta);
    }
    
    double Delta(double r) const {
        return r * r - 2.0 * M_ * r + a_ * a_;
    }

public:
    MetricaKerr(double massa_solar = 1.0, double spin = 0.0) {
        massa_kg_ = massa_solar * MASSA_SOL;
        M_ = G * massa_kg_ / C2;
        spin_ = std::max(-0.998, std::min(0.998, spin));  // Limite físico
        a_ = spin_ * M_;
    }

    // ============================================================
    // COMPONENTES DA MÉTRICA (Boyer-Lindquist)
    // ============================================================
    
    // g_tt = -(1 - 2Mr/Σ)
    double g_tt(double r, double theta) const {
        double sigma = Sigma(r, theta);
        return -(1.0 - 2.0 * M_ * r / sigma);
    }

    // g_tφ = -2Mar sin²θ / Σ
    double g_t_phi(double r, double theta) const {
        double sigma = Sigma(r, theta);
        double s = std::sin(theta);
        return -2.0 * M_ * a_ * r * s * s / sigma;
    }

    // g_rr = Σ/Δ
    double g_rr(double r, double theta) const {
        double sigma = Sigma(r, theta);
        double delta = Delta(r);
        if (std::abs(delta) < 1e-10) return 1e10;
        return sigma / delta;
    }

    // g_θθ = Σ
    double g_theta_theta(double r, double theta) const {
        return Sigma(r, theta);
    }

    // g_φφ = (r² + a² + 2Ma²r sin²θ/Σ) sin²θ
    double g_phi_phi(double r, double theta) const {
        double sigma = Sigma(r, theta);
        double s = std::sin(theta);
        return (r*r + a_*a_ + 2.0*M_*a_*a_*r*s*s/sigma) * s * s;
    }

    // ============================================================
    // HORIZONTES E ERGOSFERA
    // ============================================================

    // Horizonte externo: r+ = M + √(M² - a²)
    double horizonte_externo() const {
        return M_ + std::sqrt(M_ * M_ - a_ * a_);
    }

    // Horizonte interno: r- = M - √(M² - a²)
    double horizonte_interno() const {
        return M_ - std::sqrt(M_ * M_ - a_ * a_);
    }

    // Superfície da ergosfera: r_ergo = M + √(M² - a²cos²θ)
    double ergosfera(double theta) const {
        double c = std::cos(theta);
        return M_ + std::sqrt(M_ * M_ - a_ * a_ * c * c);
    }

    // ============================================================
    // PROPRIEDADES FÍSICAS
    // ============================================================

    // Velocidade angular do horizonte: Ω_H = a / (r+² + a²)
    double velocidade_angular_horizonte() const {
        double rp = horizonte_externo();
        return a_ / (rp * rp + a_ * a_);
    }

    // Temperatura de Hawking (Kerr)
    double temperatura_hawking() const {
        double rp = horizonte_externo();
        double rm = horizonte_interno();
        return (rp - rm) / (4.0 * M_PI * (rp * rp + a_ * a_));
    }

    // Entropia (proporcional à área)
    double entropia() const {
        double rp = horizonte_externo();
        double area = 4.0 * M_PI * (rp * rp + a_ * a_);
        return K_BOLTZMANN * C * C * C * area / (4.0 * G * H_BARRA);
    }

    // ISCO (órbita circular mais interna estável)
    double raio_isco(bool progrado = true) const {
        double z1 = 1.0 + std::cbrt(1.0 - spin_*spin_) * 
            (std::cbrt(1.0 + spin_) + std::cbrt(1.0 - spin_));
        double z2 = std::sqrt(3.0*spin_*spin_ + z1*z1);
        
        if (progrado) {
            return M_ * (3.0 + z2 - std::sqrt((3.0-z1)*(3.0+z1+2.0*z2)));
        } else {
            return M_ * (3.0 + z2 + std::sqrt((3.0-z1)*(3.0+z1+2.0*z2)));
        }
    }

    // Frame dragging angular velocity
    double frame_dragging(double r, double theta) const {
        double sigma = Sigma(r, theta);
        double A = (r*r + a_*a_)*(r*r + a_*a_) - Delta(r)*a_*a_*std::sin(theta)*std::sin(theta);
        return 2.0 * M_ * a_ * r / A;
    }

    // Getters
    double massa_geometrica() const { return M_; }
    double parametro_spin() const { return a_; }
    double spin_adimensional() const { return spin_; }
    double massa_kg() const { return massa_kg_; }
};

} // namespace BuracoNegro

#endif // KERR_HPP
