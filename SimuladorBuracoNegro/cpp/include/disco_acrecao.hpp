// ============================================================
// SimuladorBuracoNegro - Disco de Acreção
// Autor: Luiz Tiago Wilcke
// 
// Modelo físico de disco de acreção Shakura-Sunyaev
// com emissão térmica e efeitos relativísticos
// ============================================================

#ifndef DISCO_ACRECAO_HPP
#define DISCO_ACRECAO_HPP

#include "constantes.hpp"
#include "schwarzschild.hpp"
#include <cmath>
#include <array>

namespace BuracoNegro {

// ============================================================
// ESTRUTURAS DE DADOS
// ============================================================

// Propriedades de um anel do disco
struct AnelDisco {
    double raio;            // Raio do anel (m)
    double temperatura;     // Temperatura local (K)
    double luminosidade;    // Luminosidade por unidade de área (W/m²)
    double velocidade_orbital;  // Velocidade Kepleriana
    double fator_redshift;  // Redshift gravitacional
};

// Cor RGB para renderização
struct CorRGB {
    double r, g, b;  // Valores 0.0 - 1.0
    
    CorRGB operator*(double s) const {
        return {r * s, g * s, b * s};
    }
    
    CorRGB operator+(const CorRGB& c) const {
        return {r + c.r, g + c.g, b + c.b};
    }
    
    void clamp() {
        r = std::max(0.0, std::min(1.0, r));
        g = std::max(0.0, std::min(1.0, g));
        b = std::max(0.0, std::min(1.0, b));
    }
};

// ============================================================
// CLASSE DISCO DE ACREÇÃO
// ============================================================

class DiscoAcrecao {
private:
    double massa_bh_;           // Massa do buraco negro (kg)
    double taxa_acrecao_;       // Taxa de acreção (kg/s)
    double raio_interno_;       // ISCO (m)
    double raio_externo_;       // Borda externa do disco (m)
    double spin_;               // Parâmetro de spin (0-1)
    
    // Constante de Stefan-Boltzmann para cálculo de temperatura
    static constexpr double SIGMA_SB = 5.670374419e-8;
    
public:
    // Construtor
    DiscoAcrecao(double massa_solar = 10.0, 
                 double taxa_acrecao_eddington = 0.1,
                 double spin = 0.0)
    {
        massa_bh_ = massa_solar * MASSA_SOL;
        spin_ = std::max(0.0, std::min(0.998, spin));
        
        // Raio de Schwarzschild
        double rs = 2.0 * G * massa_bh_ / C2;
        
        // ISCO depende do spin
        if (spin_ < 0.01) {
            raio_interno_ = 3.0 * rs;  // Schwarzschild
        } else {
            // Aproximação para Kerr
            raio_interno_ = rs * (3.0 + spin_ - std::sqrt((3.0 - spin_) * (1.0 + spin_)));
        }
        
        // Raio externo típico = 1000 rs
        raio_externo_ = 500.0 * rs;
        
        // Taxa de acreção em termos de Eddington
        // L_Edd = 4πGMm_p c / σ_T ≈ 1.26×10³⁸ M/M_sol W
        double luminosidade_eddington = 1.26e38 * massa_solar;  // W
        double eficiencia = 0.1;  // ~10% para Schwarzschild
        taxa_acrecao_ = taxa_acrecao_eddington * luminosidade_eddington / (eficiencia * C2);
    }

    // ============================================================
    // PERFIL DE TEMPERATURA SHAKURA-SUNYAEV
    // ============================================================
    
    // Temperatura efetiva em função do raio
    // T(r) = T_* × (r/r_in)^(-3/4) × f(r)
    // onde f(r) = [1 - (r_in/r)^(1/2)]^(1/4)
    double temperatura(double raio) const {
        if (raio < raio_interno_ || raio > raio_externo_) {
            return 0.0;
        }
        
        // Temperatura característica no ISCO
        // T_* = [3GMṀ / (8πσr³)]^(1/4)
        double T_estrela = std::pow(
            3.0 * G * massa_bh_ * taxa_acrecao_ / 
            (8.0 * M_PI * SIGMA_SB * std::pow(raio_interno_, 3)),
            0.25
        );
        
        // Perfil radial
        double x = raio / raio_interno_;
        double fator_radial = std::pow(x, -0.75);
        double fator_borda = std::pow(1.0 - std::sqrt(1.0 / x), 0.25);
        
        return T_estrela * fator_radial * fator_borda;
    }

    // ============================================================
    // EMISSÃO DE CORPO NEGRO
    // ============================================================
    
    // Lei de Planck: B(ν,T) = (2hν³/c²) × 1/(exp(hν/kT) - 1)
    double planck(double frequencia, double T) const {
        if (T <= 0.0) return 0.0;
        
        double x = H_PLANCK * frequencia / (K_BOLTZMANN * T);
        if (x > 700.0) return 0.0;  // Evita overflow
        
        return (2.0 * H_PLANCK * std::pow(frequencia, 3) / C2) / 
               (std::exp(x) - 1.0);
    }
    
    // Cor aproximada de corpo negro (RGB)
    CorRGB cor_corpo_negro(double T) const {
        // Algoritmo baseado na aproximação de Planck
        // para conversão temperatura → cor visível
        
        if (T <= 0.0) return {0.0, 0.0, 0.0};
        
        double r, g, b;
        
        // Normaliza temperatura para escala prática
        double t = T / 100.0;
        
        // Canal vermelho
        if (t <= 66.0) {
            r = 1.0;
        } else {
            r = 1.29293618606274 * std::pow(t - 60.0, -0.1332047592);
        }
        
        // Canal verde
        if (t <= 66.0) {
            g = 0.390081578769871 * std::log(t) - 0.631841443788627;
        } else {
            g = 1.12989086089529 * std::pow(t - 60.0, -0.0755148492);
        }
        
        // Canal azul
        if (t >= 66.0) {
            b = 1.0;
        } else if (t <= 19.0) {
            b = 0.0;
        } else {
            b = 0.543206789110196 * std::log(t - 10.0) - 1.19625408914;
        }
        
        CorRGB cor = {r, g, b};
        cor.clamp();
        return cor;
    }

    // ============================================================
    // EFEITOS RELATIVÍSTICOS
    // ============================================================
    
    // Fator de redshift gravitacional
    double fator_redshift(double raio) const {
        double rs = 2.0 * G * massa_bh_ / C2;
        if (raio <= rs) return 0.0;
        return std::sqrt(1.0 - rs / raio);
    }
    
    // Velocidade orbital Kepleriana
    double velocidade_kepleriana(double raio) const {
        return std::sqrt(G * massa_bh_ / raio);
    }
    
    // Doppler beaming relativístico
    // D = 1 / γ(1 - β·n̂)  onde β = v/c
    double fator_doppler(double raio, double angulo_observador) const {
        double v = velocidade_kepleriana(raio);
        double beta = v / C;
        double gamma = 1.0 / std::sqrt(1.0 - beta * beta);
        
        // Ângulo entre velocidade e direção do observador
        double cos_phi = std::cos(angulo_observador);
        
        return 1.0 / (gamma * (1.0 - beta * cos_phi));
    }
    
    // Intensidade observada com efeitos relativísticos
    // I_obs = D⁴ × I_emitido (para emissão isotrópica)
    CorRGB intensidade_observada(double raio, double angulo_obs) const {
        double T = temperatura(raio);
        if (T <= 0.0) return {0.0, 0.0, 0.0};
        
        CorRGB cor = cor_corpo_negro(T);
        
        // Efeitos relativísticos
        double D = fator_doppler(raio, angulo_obs);
        double z = fator_redshift(raio);
        
        // Combinação de Doppler e redshift gravitacional
        double fator_total = std::pow(D * z, 4);
        
        return cor * fator_total;
    }

    // ============================================================
    // PROPRIEDADES DO DISCO
    // ============================================================
    
    AnelDisco anel(double raio) const {
        AnelDisco a;
        a.raio = raio;
        a.temperatura = temperatura(raio);
        a.luminosidade = SIGMA_SB * std::pow(a.temperatura, 4);
        a.velocidade_orbital = velocidade_kepleriana(raio);
        a.fator_redshift = fator_redshift(raio);
        return a;
    }
    
    // Luminosidade total do disco
    // L = ηṀc² onde η ≈ 1 - √(r_isco/r) para discos finos
    double luminosidade_total() const {
        double rs = 2.0 * G * massa_bh_ / C2;
        double eta = 1.0 - std::sqrt(rs / raio_interno_);
        return eta * taxa_acrecao_ * C2;
    }
    
    // Temperatura máxima (próximo ao ISCO)
    double temperatura_maxima() const {
        // Máximo ocorre em r ≈ 49/36 × r_isco
        double r_max = 1.361 * raio_interno_;
        return temperatura(r_max);
    }
    
    // Getters
    double raio_interno() const { return raio_interno_; }
    double raio_externo() const { return raio_externo_; }
    double massa_bh() const { return massa_bh_; }
    double taxa_acrecao() const { return taxa_acrecao_; }
    
    // Verifica se ponto está no disco
    bool no_disco(double raio) const {
        return raio >= raio_interno_ && raio <= raio_externo_;
    }
};

} // namespace BuracoNegro

#endif // DISCO_ACRECAO_HPP
