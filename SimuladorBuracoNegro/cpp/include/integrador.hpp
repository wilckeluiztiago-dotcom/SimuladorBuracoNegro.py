// ============================================================
// SimuladorBuracoNegro - Integrador Geodésico
// Autor: Luiz Tiago Wilcke
// 
// Integração numérica de geodésicas usando RK4
// ============================================================

#ifndef INTEGRADOR_HPP
#define INTEGRADOR_HPP

#include "schwarzschild.hpp"
#include <vector>
#include <functional>

namespace BuracoNegro {

// Estado de integração (8 variáveis)
struct EstadoGeodesica {
    double t, r, theta, phi;       // Posição
    double u_t, u_r, u_theta, u_phi; // Quadrivelocidade
    
    std::array<double, 8> como_array() const {
        return {t, r, theta, phi, u_t, u_r, u_theta, u_phi};
    }
    
    static EstadoGeodesica de_array(const std::array<double, 8>& arr) {
        return {arr[0], arr[1], arr[2], arr[3], arr[4], arr[5], arr[6], arr[7]};
    }
};

// Ponto na trajetória para visualização
struct PontoTrajetoria {
    double lambda;  // Parâmetro afim
    double t, r, theta, phi;
    double x, y, z; // Coordenadas cartesianas
};

class IntegradorGeodesico {
private:
    MetricaSchwarzschild metrica_;
    double passo_;
    double r_min_;
    
public:
    IntegradorGeodesico(const MetricaSchwarzschild& metrica, 
                        double passo = 0.01)
        : metrica_(metrica), passo_(passo) {
        r_min_ = metrica_.raio_schwarzschild() * 1.001;
    }
    
    // Passo RK4
    EstadoGeodesica passo_rk4(const EstadoGeodesica& estado) const {
        auto arr = estado.como_array();
        
        auto derivs = [this](const std::array<double, 8>& y) {
            Coordenadas pos{y[0], y[1], y[2], y[3]};
            Quadrivelocidade vel{y[4], y[5], y[6], y[7]};
            return metrica_.derivadas_geodesica(pos, vel);
        };
        
        auto k1 = derivs(arr);
        
        std::array<double, 8> y2;
        for (int i = 0; i < 8; i++) y2[i] = arr[i] + 0.5 * passo_ * k1[i];
        auto k2 = derivs(y2);
        
        std::array<double, 8> y3;
        for (int i = 0; i < 8; i++) y3[i] = arr[i] + 0.5 * passo_ * k2[i];
        auto k3 = derivs(y3);
        
        std::array<double, 8> y4;
        for (int i = 0; i < 8; i++) y4[i] = arr[i] + passo_ * k3[i];
        auto k4 = derivs(y4);
        
        std::array<double, 8> novo;
        for (int i = 0; i < 8; i++) {
            novo[i] = arr[i] + passo_ * (k1[i] + 2*k2[i] + 2*k3[i] + k4[i]) / 6.0;
        }
        
        return EstadoGeodesica::de_array(novo);
    }
    
    // Integrar trajetória completa
    std::vector<PontoTrajetoria> integrar(
        const EstadoGeodesica& inicial,
        double lambda_max,
        int max_pontos = 10000
    ) {
        std::vector<PontoTrajetoria> trajetoria;
        EstadoGeodesica estado = inicial;
        double lambda = 0.0;
        
        while (lambda < lambda_max && trajetoria.size() < max_pontos) {
            // Verifica se caiu no horizonte
            if (estado.r < r_min_) break;
            
            // Converte para cartesianas
            double x = estado.r * std::sin(estado.theta) * std::cos(estado.phi);
            double y = estado.r * std::sin(estado.theta) * std::sin(estado.phi);
            double z = estado.r * std::cos(estado.theta);
            
            trajetoria.push_back({lambda, estado.t, estado.r, 
                                  estado.theta, estado.phi, x, y, z});
            
            estado = passo_rk4(estado);
            lambda += passo_;
        }
        
        return trajetoria;
    }
    
    // Criar condições iniciais para fóton
    EstadoGeodesica foton_inicial(
        double r, double theta, double phi,
        double direcao_r, double direcao_theta, double direcao_phi
    ) const {
        // Normaliza direção para geodésica nula
        double rs = metrica_.raio_schwarzschild();
        double f = 1.0 - rs / r;
        
        // Para geodésica nula: g_μν u^μ u^ν = 0
        double u_t = 1.0 / f;  // Energia = 1
        
        return {0.0, r, theta, phi, 
                u_t, direcao_r, direcao_theta, direcao_phi};
    }
    
    void set_passo(double p) { passo_ = p; }
    double get_passo() const { return passo_; }
};

} // namespace BuracoNegro

#endif // INTEGRADOR_HPP
