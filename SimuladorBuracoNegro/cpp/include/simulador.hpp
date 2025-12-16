// ============================================================
// SimuladorBuracoNegro - Classe Principal do Simulador
// Autor: Luiz Tiago Wilcke
// 
// Encapsula toda a simulação e gerencia configurações
// ============================================================

#ifndef SIMULADOR_HPP
#define SIMULADOR_HPP

#include "constantes.hpp"
#include "schwarzschild.hpp"
#include "kerr.hpp"
#include "disco_acrecao.hpp"
#include "ray_tracer.hpp"
#include "integrador.hpp"
#include <string>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <chrono>
#include <filesystem>

namespace BuracoNegro {

// ============================================================
// ENUMERAÇÕES
// ============================================================

enum class TipoBuracoNegro {
    SCHWARZSCHILD,
    KERR
};

enum class ModoSimulacao {
    RAY_TRACING,        // Renderização de imagem
    GEODESICAS,         // Trajetórias de partículas
    ANALISE             // Análise física
};

// ============================================================
// CONFIGURAÇÃO DA SIMULAÇÃO
// ============================================================

struct ConfiguracaoSimulacao {
    // Parâmetros do buraco negro
    TipoBuracoNegro tipo = TipoBuracoNegro::SCHWARZSCHILD;
    double massa_solar = 10.0;          // Massas solares
    double spin = 0.0;                   // Parâmetro a/M (0-0.998 para Kerr)
    
    // Parâmetros do disco de acreção
    bool incluir_disco = true;
    double taxa_eddington = 0.1;         // Fração da taxa de Eddington
    
    // Parâmetros de visualização
    int largura = 800;
    int altura = 600;
    double distancia_observador = 100.0; // Em raios de Schwarzschild
    double angulo_inclinacao = 75.0;     // Graus (0 = de cima, 90 = de lado)
    double fov = 45.0;                   // Campo de visão em graus
    
    // Parâmetros de integração
    double passo_integracao = 0.1;
    int max_passos = 10000;
    
    // Parâmetros de performance
    int num_threads = 4;
    
    // Saída
    std::string diretorio_saida = "../saida";
    std::string prefixo_arquivo = "buraco_negro";
};

// ============================================================
// CLASSE SIMULADOR
// ============================================================

class SimuladorBuracoNegro {
private:
    ConfiguracaoSimulacao config_;
    std::unique_ptr<RayTracer> ray_tracer_;
    
    // Métricas de desempenho
    double tempo_ultima_renderizacao_ = 0.0;
    
public:
    SimuladorBuracoNegro() = default;
    
    explicit SimuladorBuracoNegro(const ConfiguracaoSimulacao& config)
        : config_(config) 
    {
        inicializar();
    }
    
    // ============================================================
    // INICIALIZAÇÃO
    // ============================================================
    
    void inicializar() {
        // Cria o ray tracer
        ray_tracer_ = std::make_unique<RayTracer>(
            config_.massa_solar,
            config_.taxa_eddington
        );
        
        // Configura câmera
        Camera cam;
        cam.largura = config_.largura;
        cam.altura = config_.altura;
        cam.r_observador = config_.distancia_observador;
        cam.theta_observador = (90.0 - config_.angulo_inclinacao) * M_PI / 180.0;
        cam.fov_horizontal = config_.fov * M_PI / 180.0;
        cam.fov_vertical = config_.fov * M_PI / 180.0 * config_.altura / config_.largura;
        
        ray_tracer_->set_camera(cam);
        ray_tracer_->set_threads(config_.num_threads);
        
        // Cria diretório de saída se não existe
        std::filesystem::create_directories(config_.diretorio_saida);
    }
    
    void set_configuracao(const ConfiguracaoSimulacao& config) {
        config_ = config;
        inicializar();
    }
    
    // ============================================================
    // EXECUÇÃO
    // ============================================================
    
    bool renderizar() {
        if (!ray_tracer_) {
            std::cerr << "Erro: Ray tracer não inicializado!" << std::endl;
            return false;
        }
        
        std::cout << "\n";
        std::cout << "╔══════════════════════════════════════════════════════════════╗\n";
        std::cout << "║         SIMULADOR DE BURACO NEGRO RELATIVÍSTICO              ║\n";
        std::cout << "╠══════════════════════════════════════════════════════════════╣\n";
        std::cout << "║  Autor: Luiz Tiago Wilcke                                    ║\n";
        std::cout << "╚══════════════════════════════════════════════════════════════╝\n";
        std::cout << "\n";
        
        imprimir_parametros();
        
        std::cout << "\n[INICIANDO RENDERIZAÇÃO...]\n\n";
        
        auto inicio = std::chrono::high_resolution_clock::now();
        
        // Executa ray tracing
        auto imagem = ray_tracer_->renderizar();
        
        auto fim = std::chrono::high_resolution_clock::now();
        tempo_ultima_renderizacao_ = 
            std::chrono::duration<double>(fim - inicio).count();
        
        std::cout << "\n[RENDERIZAÇÃO COMPLETA]\n";
        std::cout << "  Tempo: " << std::fixed << std::setprecision(2) 
                  << tempo_ultima_renderizacao_ << " segundos\n";
        std::cout << "  Pixels: " << config_.largura * config_.altura << "\n";
        std::cout << "  Taxa: " << std::fixed << std::setprecision(0)
                  << (config_.largura * config_.altura) / tempo_ultima_renderizacao_
                  << " pixels/s\n\n";
        
        // Salva imagem
        std::string nome_arquivo = gerar_nome_arquivo();
        std::string caminho_ppm = config_.diretorio_saida + "/" + nome_arquivo + ".ppm";
        
        if (ray_tracer_->salvar_ppm(imagem, caminho_ppm)) {
            std::cout << "[SALVO] " << caminho_ppm << "\n";
        } else {
            std::cerr << "[ERRO] Não foi possível salvar a imagem!\n";
            return false;
        }
        
        return true;
    }
    
    // ============================================================
    // ANÁLISE FÍSICA
    // ============================================================
    
    void imprimir_parametros() const {
        double rs = raio_schwarzschild_solar(config_.massa_solar);
        
        std::cout << "┌─────────────────────────────────────────────────────────────┐\n";
        std::cout << "│ PARÂMETROS DO BURACO NEGRO                                  │\n";
        std::cout << "├─────────────────────────────────────────────────────────────┤\n";
        
        std::cout << "│ Tipo: ";
        if (config_.tipo == TipoBuracoNegro::SCHWARZSCHILD) {
            std::cout << "Schwarzschild (estático, esférico)                     │\n";
        } else {
            std::cout << "Kerr (rotativo), spin = " << std::fixed 
                      << std::setprecision(3) << config_.spin << "                   │\n";
        }
        
        std::cout << std::fixed << std::setprecision(1);
        std::cout << "│ Massa: " << config_.massa_solar << " M☉ ("
                  << std::scientific << std::setprecision(3) 
                  << config_.massa_solar * MASSA_SOL << " kg)               │\n";
        
        std::cout << std::fixed << std::setprecision(3);
        std::cout << "│ Raio de Schwarzschild: " << rs << " m ("
                  << std::fixed << std::setprecision(1) 
                  << rs / 1000.0 << " km)                    │\n";
        
        double T_hawking = temperatura_hawking(config_.massa_solar * MASSA_SOL);
        std::cout << std::scientific << std::setprecision(2);
        std::cout << "│ Temperatura de Hawking: " << T_hawking << " K                        │\n";
        
        std::cout << "├─────────────────────────────────────────────────────────────┤\n";
        std::cout << "│ CONFIGURAÇÃO DA CÂMERA                                      │\n";
        std::cout << "├─────────────────────────────────────────────────────────────┤\n";
        std::cout << std::fixed << std::setprecision(0);
        std::cout << "│ Resolução: " << config_.largura << " x " << config_.altura 
                  << " pixels                                     │\n";
        std::cout << "│ Distância: " << config_.distancia_observador << " rs ("
                  << std::scientific << std::setprecision(2)
                  << config_.distancia_observador * rs << " m)            │\n";
        std::cout << std::fixed << std::setprecision(1);
        std::cout << "│ Inclinação: " << config_.angulo_inclinacao << "°                                               │\n";
        std::cout << "│ Campo de visão: " << config_.fov << "°                                           │\n";
        std::cout << "│ Threads: " << config_.num_threads << "                                                     │\n";
        std::cout << "└─────────────────────────────────────────────────────────────┘\n";
    }
    
    void analise_fisica() const {
        double M = config_.massa_solar * MASSA_SOL;
        double rs = raio_schwarzschild(M);
        
        std::cout << "\n";
        std::cout << "╔══════════════════════════════════════════════════════════════╗\n";
        std::cout << "║              ANÁLISE FÍSICA DO BURACO NEGRO                  ║\n";
        std::cout << "╚══════════════════════════════════════════════════════════════╝\n\n";
        
        std::cout << "PROPRIEDADES GEOMÉTRICAS:\n";
        std::cout << "  • Raio de Schwarzschild: " << std::scientific << rs << " m\n";
        std::cout << "  • Raio ISCO: " << 3.0 * rs << " m\n";
        std::cout << "  • Raio da esfera de fótons: " << 1.5 * rs << " m\n";
        std::cout << "  • Área do horizonte: " << 4.0 * M_PI * rs * rs << " m²\n\n";
        
        std::cout << "PROPRIEDADES TERMODINÂMICAS:\n";
        double T = temperatura_hawking(M);
        double S = entropia_bekenstein_hawking(M);
        double L = luminosidade_hawking(M);
        double t_evap = tempo_evaporacao(M);
        
        std::cout << "  • Temperatura de Hawking: " << T << " K\n";
        std::cout << "  • Entropia (Bekenstein-Hawking): " << S << " J/K\n";
        std::cout << "  • Luminosidade de Hawking: " << L << " W\n";
        std::cout << "  • Tempo de evaporação: " << t_evap << " s ("
                  << t_evap / (365.25 * 24 * 3600) << " anos)\n\n";
        
        std::cout << "EFEITOS RELATIVÍSTICOS:\n";
        for (double r_fator : {1.5, 2.0, 3.0, 5.0, 10.0, 100.0}) {
            double r = r_fator * rs;
            double dilatacao = std::sqrt(1.0 - rs / r);
            double v_escape = C * std::sqrt(rs / r);
            
            std::cout << "  r = " << std::fixed << std::setprecision(1) << r_fator 
                      << " rs: dilatação temporal = " << std::setprecision(4) << dilatacao
                      << ", v_escape = " << std::scientific << std::setprecision(2) 
                      << v_escape << " m/s (" 
                      << std::fixed << std::setprecision(1) << 100.0 * v_escape / C << "% c)\n";
        }
    }
    
    // ============================================================
    // UTILITÁRIOS
    // ============================================================
    
    std::string gerar_nome_arquivo() const {
        std::ostringstream ss;
        ss << config_.prefixo_arquivo;
        ss << "_M" << static_cast<int>(config_.massa_solar);
        ss << "_inc" << static_cast<int>(config_.angulo_inclinacao);
        ss << "_" << config_.largura << "x" << config_.altura;
        
        // Adiciona timestamp
        auto now = std::chrono::system_clock::now();
        auto time_t = std::chrono::system_clock::to_time_t(now);
        ss << "_" << std::put_time(std::localtime(&time_t), "%Y%m%d_%H%M%S");
        
        return ss.str();
    }
    
    // Getters
    const ConfiguracaoSimulacao& configuracao() const { return config_; }
    double tempo_renderizacao() const { return tempo_ultima_renderizacao_; }
};

} // namespace BuracoNegro

#endif // SIMULADOR_HPP
