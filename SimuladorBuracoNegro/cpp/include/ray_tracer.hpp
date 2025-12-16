// ============================================================
// SimuladorBuracoNegro - Ray Tracer Relativístico
// Autor: Luiz Tiago Wilcke
// 
// Traçado de raios em espaço-tempo curvo
// para visualização de buracos negros
// ============================================================

#ifndef RAY_TRACER_HPP
#define RAY_TRACER_HPP

#include "constantes.hpp"
#include "schwarzschild.hpp"
#include "kerr.hpp"
#include "disco_acrecao.hpp"
#include "integrador.hpp"
#include <cmath>
#include <array>
#include <vector>
#include <fstream>
#include <iomanip>
#include <thread>
#include <mutex>
#include <atomic>

namespace BuracoNegro {

// ============================================================
// ESTRUTURAS DE DADOS
// ============================================================

// Pixel da imagem
struct Pixel {
    double r, g, b;
    
    Pixel() : r(0), g(0), b(0) {}
    Pixel(double _r, double _g, double _b) : r(_r), g(_g), b(_b) {}
    
    Pixel operator+(const Pixel& p) const {
        return Pixel(r + p.r, g + p.g, b + p.b);
    }
    
    Pixel operator*(double s) const {
        return Pixel(r * s, g * s, b * s);
    }
    
    void clamp() {
        r = std::max(0.0, std::min(1.0, r));
        g = std::max(0.0, std::min(1.0, g));
        b = std::max(0.0, std::min(1.0, b));
    }
    
    uint8_t r_byte() const { return static_cast<uint8_t>(r * 255); }
    uint8_t g_byte() const { return static_cast<uint8_t>(g * 255); }
    uint8_t b_byte() const { return static_cast<uint8_t>(b * 255); }
};

// Configuração da câmera
struct Camera {
    double r_observador;    // Distância do observador ao BH
    double theta_observador; // Ângulo polar (0 = face-on, π/2 = edge-on)
    double fov_horizontal;  // Campo de visão horizontal (radianos)
    double fov_vertical;    // Campo de visão vertical (radianos)
    int largura;            // Resolução horizontal
    int altura;             // Resolução vertical
    
    Camera(double r = 100.0, double theta = M_PI/3.0, 
           double fov_h = M_PI/4.0, double fov_v = M_PI/4.0,
           int w = 800, int h = 600)
        : r_observador(r), theta_observador(theta),
          fov_horizontal(fov_h), fov_vertical(fov_v),
          largura(w), altura(h) {}
};

// Resultado de um raio traçado
struct ResultadoRaio {
    enum Destino { HORIZONTE, DISCO, INFINITO, ERRO };
    Destino destino;
    double r_impacto;       // Raio onde atingiu (disco) ou parou
    double phi_impacto;     // Ângulo phi no impacto
    double theta_impacto;   // Ângulo theta no impacto
    int passos;             // Número de passos de integração
    Pixel cor;              // Cor resultante
};

// ============================================================
// CLASSE RAY TRACER
// ============================================================

class RayTracer {
private:
    MetricaSchwarzschild metrica_;
    DiscoAcrecao disco_;
    Camera camera_;
    double rs_;             // Raio de Schwarzschild
    
    // Parâmetros de integração
    double passo_inicial_ = 0.1;
    int max_passos_ = 10000;
    double tolerancia_horizonte_ = 1.001;
    
    // Textura de fundo (grid celestial)
    bool usar_fundo_ = true;
    
    // Multithreading
    int num_threads_ = 4;
    std::mutex mutex_progresso_;
    std::atomic<int> linhas_processadas_{0};

public:
    RayTracer(double massa_solar = 10.0, double taxa_eddington = 0.1)
        : metrica_(massa_solar),
          disco_(massa_solar, taxa_eddington),
          camera_()
    {
        rs_ = metrica_.raio_schwarzschild();
        camera_.r_observador = 100.0 * rs_;
    }

    // ============================================================
    // CONFIGURAÇÃO
    // ============================================================
    
    void set_camera(const Camera& cam) {
        camera_ = cam;
        camera_.r_observador *= rs_;  // Converte para unidades de rs
    }
    
    void set_resolucao(int largura, int altura) {
        camera_.largura = largura;
        camera_.altura = altura;
    }
    
    void set_angulo_observador(double theta) {
        camera_.theta_observador = theta;
    }
    
    void set_threads(int n) {
        num_threads_ = std::max(1, n);
    }

    // ============================================================
    // TRAÇADO DE RAIO INDIVIDUAL
    // ============================================================
    
    ResultadoRaio tracar_raio(double alfa, double beta) const {
        ResultadoRaio resultado;
        resultado.destino = ResultadoRaio::INFINITO;
        resultado.passos = 0;
        
        // Posição inicial do observador
        double r0 = camera_.r_observador;
        double theta0 = camera_.theta_observador;
        double phi0 = 0.0;
        double t0 = 0.0;
        
        // Direção inicial do fóton (reversa, do observador para o BH)
        // alfa = desvio horizontal, beta = desvio vertical
        double f = 1.0 - rs_ / r0;
        
        // Componentes da quadrivelocidade para fóton
        double u_t = 1.0 / f;  // Energia conservada
        double u_theta = beta / r0;
        double u_phi = alfa / (r0 * std::sin(theta0));
        
        // u_r calculado pela condição de geodésica nula: g_μν u^μ u^ν = 0
        // Para Schwarzschild: -f(u_t)² + (1/f)(u_r)² + r²(u_θ)² + r²sin²θ(u_φ)² = 0
        double termo1 = f * u_t * u_t;
        double termo2 = r0 * r0 * u_theta * u_theta;
        double termo3 = r0 * r0 * std::sin(theta0) * std::sin(theta0) * u_phi * u_phi;
        double u_r_sq = f * (termo1 - termo2 - termo3);
        
        // Raio para dentro (negativo)
        double u_r = -std::sqrt(std::max(0.0, u_r_sq));
        
        // Estado inicial
        EstadoGeodesica estado = {t0, r0, theta0, phi0, u_t, u_r, u_theta, u_phi};
        
        // Integrador
        IntegradorGeodesico integrador(metrica_, passo_inicial_);
        
        // Loop de integração
        double r_anterior = r0;
        for (int i = 0; i < max_passos_; i++) {
            resultado.passos = i;
            
            // Verifica se atingiu o horizonte
            if (estado.r < rs_ * tolerancia_horizonte_) {
                resultado.destino = ResultadoRaio::HORIZONTE;
                resultado.r_impacto = estado.r;
                resultado.cor = Pixel(0.0, 0.0, 0.0);  // Negro absoluto
                return resultado;
            }
            
            // Verifica se passou pelo plano do disco (θ ≈ π/2)
            double theta_disco = M_PI / 2.0;
            double delta_theta = std::abs(estado.theta - theta_disco);
            
            if (delta_theta < 0.01 && disco_.no_disco(estado.r)) {
                resultado.destino = ResultadoRaio::DISCO;
                resultado.r_impacto = estado.r;
                resultado.phi_impacto = estado.phi;
                resultado.theta_impacto = estado.theta;
                
                // Cor do disco
                CorRGB cor_disco = disco_.intensidade_observada(estado.r, estado.phi);
                resultado.cor = Pixel(cor_disco.r, cor_disco.g, cor_disco.b);
                return resultado;
            }
            
            // Verifica se escapou para infinito
            if (estado.r > camera_.r_observador * 2.0) {
                resultado.destino = ResultadoRaio::INFINITO;
                
                if (usar_fundo_) {
                    // Grid celestial de fundo
                    resultado.cor = cor_fundo(estado.theta, estado.phi);
                } else {
                    resultado.cor = Pixel(0.02, 0.02, 0.05);  // Azul escuro
                }
                return resultado;
            }
            
            // Passo adaptativo baseado na distância
            double fator_passo = std::sqrt(estado.r / rs_);
            integrador.set_passo(passo_inicial_ * fator_passo);
            
            // Integra um passo
            r_anterior = estado.r;
            estado = integrador.passo_rk4(estado);
            
            // Mantém θ em [0, π]
            if (estado.theta < 0.0) {
                estado.theta = -estado.theta;
                estado.u_theta = -estado.u_theta;
            }
            if (estado.theta > M_PI) {
                estado.theta = 2.0 * M_PI - estado.theta;
                estado.u_theta = -estado.u_theta;
            }
        }
        
        resultado.destino = ResultadoRaio::ERRO;
        resultado.cor = Pixel(1.0, 0.0, 1.0);  // Magenta = erro
        return resultado;
    }

    // ============================================================
    // COR DE FUNDO (GRID CELESTIAL)
    // ============================================================
    
    Pixel cor_fundo(double theta, double phi) const {
        // Normaliza ângulos
        phi = std::fmod(phi, 2.0 * M_PI);
        if (phi < 0) phi += 2.0 * M_PI;
        
        // Grid com linhas a cada 15°
        double espessura = 0.02;
        double lat = theta - M_PI / 2.0;  // Latitude: -π/2 a π/2
        double lon = phi;                  // Longitude: 0 a 2π
        
        // Linhas de latitude (constante θ)
        bool linha_lat = false;
        for (double l = -M_PI / 2.0; l <= M_PI / 2.0; l += M_PI / 12.0) {
            if (std::abs(lat - l) < espessura) {
                linha_lat = true;
                break;
            }
        }
        
        // Linhas de longitude (constante φ)
        bool linha_lon = false;
        for (double l = 0.0; l < 2.0 * M_PI; l += M_PI / 12.0) {
            double diff = std::abs(lon - l);
            if (diff < espessura || (2.0 * M_PI - diff) < espessura) {
                linha_lon = true;
                break;
            }
        }
        
        // Cor das linhas
        if (linha_lat || linha_lon) {
            // Gradiente azul-roxo baseado na posição
            double h = (lon / (2.0 * M_PI));  // Hue 0-1
            return Pixel(0.2 + 0.3 * h, 0.1, 0.4 + 0.2 * (1.0 - h));
        }
        
        // Fundo com estrelas (ruído pseudo-aleatório)
        double seed = theta * 100.0 + phi * 57.0;
        double estrela = std::sin(seed * 12345.6789);
        estrela = (estrela + 1.0) / 2.0;
        estrela = std::pow(estrela, 100.0);  // Poucas estrelas brilhantes
        
        return Pixel(0.01 + 0.5 * estrela, 0.01 + 0.5 * estrela, 0.03 + 0.5 * estrela);
    }

    // ============================================================
    // RENDERIZAÇÃO COMPLETA
    // ============================================================
    
    std::vector<std::vector<Pixel>> renderizar() {
        int largura = camera_.largura;
        int altura = camera_.altura;
        
        std::vector<std::vector<Pixel>> imagem(altura, std::vector<Pixel>(largura));
        linhas_processadas_ = 0;
        
        // Função para processar um bloco de linhas
        auto processar_bloco = [&](int linha_inicio, int linha_fim) {
            for (int j = linha_inicio; j < linha_fim; j++) {
                for (int i = 0; i < largura; i++) {
                    // Converte pixel para coordenadas de impacto
                    double alfa = (i - largura / 2.0) / largura * camera_.fov_horizontal * camera_.r_observador;
                    double beta = (j - altura / 2.0) / altura * camera_.fov_vertical * camera_.r_observador;
                    
                    ResultadoRaio resultado = tracar_raio(alfa, beta);
                    imagem[j][i] = resultado.cor;
                }
                
                linhas_processadas_++;
            }
        };
        
        // Divide trabalho entre threads
        std::vector<std::thread> threads;
        int linhas_por_thread = altura / num_threads_;
        
        for (int t = 0; t < num_threads_; t++) {
            int inicio = t * linhas_por_thread;
            int fim = (t == num_threads_ - 1) ? altura : (t + 1) * linhas_por_thread;
            threads.emplace_back(processar_bloco, inicio, fim);
        }
        
        // Aguarda todas as threads
        for (auto& th : threads) {
            th.join();
        }
        
        return imagem;
    }
    
    // Progresso da renderização (0.0 - 1.0)
    double progresso() const {
        return static_cast<double>(linhas_processadas_) / camera_.altura;
    }

    // ============================================================
    // EXPORTAÇÃO DE IMAGEM
    // ============================================================
    
    // Salva como PPM (formato simples)
    bool salvar_ppm(const std::vector<std::vector<Pixel>>& imagem, 
                    const std::string& caminho) const {
        std::ofstream arquivo(caminho, std::ios::binary);
        if (!arquivo.is_open()) return false;
        
        int altura = imagem.size();
        int largura = imagem[0].size();
        
        arquivo << "P6\n" << largura << " " << altura << "\n255\n";
        
        for (int j = 0; j < altura; j++) {
            for (int i = 0; i < largura; i++) {
                Pixel p = imagem[j][i];
                p.clamp();
                arquivo << p.r_byte() << p.g_byte() << p.b_byte();
            }
        }
        
        arquivo.close();
        return true;
    }
    
    // Salva como CSV (para análise)
    bool salvar_csv(const std::vector<std::vector<Pixel>>& imagem,
                    const std::string& caminho) const {
        std::ofstream arquivo(caminho);
        if (!arquivo.is_open()) return false;
        
        arquivo << "x,y,r,g,b\n";
        
        int altura = imagem.size();
        int largura = imagem[0].size();
        
        for (int j = 0; j < altura; j++) {
            for (int i = 0; i < largura; i++) {
                const Pixel& p = imagem[j][i];
                arquivo << i << "," << j << ","
                       << std::fixed << std::setprecision(4)
                       << p.r << "," << p.g << "," << p.b << "\n";
            }
        }
        
        arquivo.close();
        return true;
    }

    // ============================================================
    // GETTERS
    // ============================================================
    
    double raio_schwarzschild() const { return rs_; }
    const Camera& camera() const { return camera_; }
    const DiscoAcrecao& disco() const { return disco_; }
};

} // namespace BuracoNegro

#endif // RAY_TRACER_HPP
