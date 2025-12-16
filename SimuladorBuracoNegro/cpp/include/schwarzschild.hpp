// ============================================================
// SimuladorBuracoNegro - Métrica de Schwarzschild
// Autor: Luiz Tiago Wilcke
// 
// Implementa a geometria do buraco negro estático esférico
// ============================================================

#ifndef SCHWARZSCHILD_HPP
#define SCHWARZSCHILD_HPP

#include "constantes.hpp"
#include <cmath>
#include <array>
#include <vector>

namespace BuracoNegro {

// ============================================================
// ESTRUTURAS DE DADOS
// ============================================================

// Coordenadas de Schwarzschild (t, r, θ, φ)
struct Coordenadas {
    double t;      // Tempo coordenado
    double r;      // Raio
    double theta;  // Ângulo polar
    double phi;    // Ângulo azimutal
};

// Quadrivelocidade (dt/dτ, dr/dτ, dθ/dτ, dφ/dτ)
struct Quadrivelocidade {
    double u_t;
    double u_r;
    double u_theta;
    double u_phi;
};

// Estado completo de uma partícula
struct EstadoParticula {
    Coordenadas posicao;
    Quadrivelocidade velocidade;
    double massa;           // Massa da partícula
    double tempo_proprio;   // Tempo próprio τ
    bool massiva;           // true para partículas, false para fótons
};

// ============================================================
// CLASSE SCHWARZSCHILD
// ============================================================

class MetricaSchwarzschild {
private:
    double M_;          // Massa do buraco negro (em unidades geométricas G=c=1)
    double rs_;         // Raio de Schwarzschild
    double massa_kg_;   // Massa em kg

public:
    // Construtor
    explicit MetricaSchwarzschild(double massa_solar = 1.0) {
        massa_kg_ = massa_solar * MASSA_SOL;
        // Em unidades geométricas: M = GM/c²
        M_ = G * massa_kg_ / C2;
        rs_ = 2.0 * M_;
    }

    // ============================================================
    // COMPONENTES DA MÉTRICA
    // ============================================================
    
    // g_tt = -(1 - rs/r)
    double g_tt(double r) const {
        if (r <= rs_) return 0.0;
        return -(1.0 - rs_ / r);
    }

    // g_rr = 1/(1 - rs/r)
    double g_rr(double r) const {
        if (r <= rs_) return 1e10;  // Singularidade
        return 1.0 / (1.0 - rs_ / r);
    }

    // g_θθ = r²
    double g_theta_theta(double r) const {
        return r * r;
    }

    // g_φφ = r²sin²θ
    double g_phi_phi(double r, double theta) const {
        return r * r * std::sin(theta) * std::sin(theta);
    }

    // ============================================================
    // SÍMBOLOS DE CHRISTOFFEL (NÃO-NULOS)
    // ============================================================

    // Γ^t_tr = Γ^t_rt = rs / (2r(r-rs))
    double christoffel_t_tr(double r) const {
        if (r <= rs_) return 0.0;
        return rs_ / (2.0 * r * (r - rs_));
    }

    // Γ^r_tt = rs(r-rs) / (2r³)
    double christoffel_r_tt(double r) const {
        if (r <= rs_) return 0.0;
        return rs_ * (r - rs_) / (2.0 * r * r * r);
    }

    // Γ^r_rr = -rs / (2r(r-rs))
    double christoffel_r_rr(double r) const {
        if (r <= rs_) return 0.0;
        return -rs_ / (2.0 * r * (r - rs_));
    }

    // Γ^r_θθ = -(r - rs)
    double christoffel_r_theta_theta(double r) const {
        return -(r - rs_);
    }

    // Γ^r_φφ = -(r - rs)sin²θ
    double christoffel_r_phi_phi(double r, double theta) const {
        double s = std::sin(theta);
        return -(r - rs_) * s * s;
    }

    // Γ^θ_rθ = Γ^θ_θr = 1/r
    double christoffel_theta_r_theta(double r) const {
        return 1.0 / r;
    }

    // Γ^θ_φφ = -sinθcosθ
    double christoffel_theta_phi_phi(double theta) const {
        return -std::sin(theta) * std::cos(theta);
    }

    // Γ^φ_rφ = Γ^φ_φr = 1/r
    double christoffel_phi_r_phi(double r) const {
        return 1.0 / r;
    }

    // Γ^φ_θφ = Γ^φ_φθ = cotθ
    double christoffel_phi_theta_phi(double theta) const {
        return 1.0 / std::tan(theta);
    }

    // ============================================================
    // EQUAÇÕES GEODÉSICAS
    // ============================================================

    // Derivadas para integração: d²x^μ/dλ² = -Γ^μ_αβ (dx^α/dλ)(dx^β/dλ)
    std::array<double, 8> derivadas_geodesica(
        const Coordenadas& pos,
        const Quadrivelocidade& vel
    ) const {
        double r = pos.r;
        double theta = pos.theta;
        
        std::array<double, 8> derivs;
        
        // Derivadas das coordenadas (velocidades)
        derivs[0] = vel.u_t;      // dt/dλ
        derivs[1] = vel.u_r;      // dr/dλ
        derivs[2] = vel.u_theta;  // dθ/dλ
        derivs[3] = vel.u_phi;    // dφ/dλ
        
        // Derivadas das velocidades (acelerações)
        // d²t/dλ²
        derivs[4] = -2.0 * christoffel_t_tr(r) * vel.u_t * vel.u_r;
        
        // d²r/dλ²
        derivs[5] = -christoffel_r_tt(r) * vel.u_t * vel.u_t
                   - christoffel_r_rr(r) * vel.u_r * vel.u_r
                   - christoffel_r_theta_theta(r) * vel.u_theta * vel.u_theta
                   - christoffel_r_phi_phi(r, theta) * vel.u_phi * vel.u_phi;
        
        // d²θ/dλ²
        derivs[6] = -2.0 * christoffel_theta_r_theta(r) * vel.u_r * vel.u_theta
                   - christoffel_theta_phi_phi(theta) * vel.u_phi * vel.u_phi;
        
        // d²φ/dλ²
        derivs[7] = -2.0 * christoffel_phi_r_phi(r) * vel.u_r * vel.u_phi
                   - 2.0 * christoffel_phi_theta_phi(theta) * vel.u_theta * vel.u_phi;
        
        return derivs;
    }

    // ============================================================
    // INTEGRAIS DE MOVIMENTO
    // ============================================================

    // Energia específica: E/m = (1 - rs/r) dt/dτ
    double energia_especifica(double r, double u_t) const {
        return (1.0 - rs_ / r) * u_t;
    }

    // Momento angular específico: L/m = r²sin²θ dφ/dτ
    double momento_angular(double r, double theta, double u_phi) const {
        double s = std::sin(theta);
        return r * r * s * s * u_phi;
    }

    // Potencial efetivo para órbitas equatoriais
    double potencial_efetivo(double r, double L, bool massiva = true) const {
        double L2 = L * L;
        if (massiva) {
            // V²_eff = (1 - rs/r)(1 + L²/r²)
            return (1.0 - rs_ / r) * (1.0 + L2 / (r * r));
        } else {
            // Para fótons: V²_eff = (1 - rs/r) L²/r²
            return (1.0 - rs_ / r) * L2 / (r * r);
        }
    }

    // ============================================================
    // PROPRIEDADES FÍSICAS
    // ============================================================

    double raio_schwarzschild() const { return rs_; }
    double massa_geometrica() const { return M_; }
    double massa_kg() const { return massa_kg_; }
    
    double raio_isco() const { return 3.0 * rs_; }
    double raio_esfera_fotons() const { return 1.5 * rs_; }
    
    double temperatura_hawking() const {
        return BuracoNegro::temperatura_hawking(massa_kg_);
    }
    
    double entropia() const {
        return entropia_bekenstein_hawking(massa_kg_);
    }

    // Fator de dilatação temporal
    double dilatacao_temporal(double r) const {
        if (r <= rs_) return 0.0;
        return std::sqrt(1.0 - rs_ / r);
    }

    // Redshift gravitacional
    double redshift(double r_emissor, double r_observador) const {
        return dilatacao_temporal(r_observador) / dilatacao_temporal(r_emissor) - 1.0;
    }

    // Velocidade de escape
    double velocidade_escape(double r) const {
        if (r <= rs_) return C;
        return C * std::sqrt(rs_ / r);
    }

    // ============================================================
    // TENSOR DE CURVATURA
    // ============================================================

    // Escalar de Kretschmann: K = R_μνρσ R^μνρσ = 48M²/r⁶
    double kretschmann(double r) const {
        return 48.0 * M_ * M_ / std::pow(r, 6);
    }

    // Tensor de Ricci (é zero para solução de vácuo)
    double ricci_scalar() const { return 0.0; }
};

} // namespace BuracoNegro

#endif // SCHWARZSCHILD_HPP
