// ============================================================
// SimuladorBuracoNegro - Programa Principal
// Autor: Luiz Tiago Wilcke
// 
// Simulador de buracos negros com ray tracing relativístico
// ============================================================

#include "simulador.hpp"
#include <cstdlib>
#include <cstring>
#include <getopt.h>

using namespace BuracoNegro;

// ============================================================
// INTERFACE DE LINHA DE COMANDO
// ============================================================

void imprimir_ajuda(const char* nome_programa) {
    std::cout << "\n";
    std::cout << "Uso: " << nome_programa << " [opções]\n\n";
    std::cout << "OPÇÕES:\n";
    std::cout << "  -m, --massa <M>        Massa em massas solares (padrão: 10)\n";
    std::cout << "  -s, --spin <a>         Parâmetro de spin 0-0.998 (padrão: 0)\n";
    std::cout << "  -i, --inclinacao <θ>   Ângulo de inclinação em graus (padrão: 75)\n";
    std::cout << "  -d, --distancia <r>    Distância do observador em rs (padrão: 100)\n";
    std::cout << "  -W, --largura <px>     Largura da imagem (padrão: 800)\n";
    std::cout << "  -H, --altura <px>      Altura da imagem (padrão: 600)\n";
    std::cout << "  -f, --fov <graus>      Campo de visão (padrão: 45)\n";
    std::cout << "  -t, --threads <n>      Número de threads (padrão: 4)\n";
    std::cout << "  -o, --saida <dir>      Diretório de saída (padrão: ../saida)\n";
    std::cout << "  -a, --analise          Apenas imprimir análise física\n";
    std::cout << "  -h, --ajuda            Mostrar esta mensagem\n";
    std::cout << "\n";
    std::cout << "EXEMPLOS:\n";
    std::cout << "  " << nome_programa << " -m 20 -i 60 -W 1920 -H 1080\n";
    std::cout << "  " << nome_programa << " -m 10 -s 0.9 -i 85 -t 8\n";
    std::cout << "  " << nome_programa << " -a -m 100  # Apenas análise de BH de 100 M☉\n";
    std::cout << "\n";
}

void menu_interativo(ConfiguracaoSimulacao& config) {
    std::cout << "\n";
    std::cout << "╔══════════════════════════════════════════════════════════════╗\n";
    std::cout << "║         SIMULADOR DE BURACO NEGRO RELATIVÍSTICO              ║\n";
    std::cout << "║                    Luiz Tiago Wilcke                         ║\n";
    std::cout << "╠══════════════════════════════════════════════════════════════╣\n";
    std::cout << "║  Este programa simula a aparência visual de um buraco negro  ║\n";
    std::cout << "║  usando ray tracing em espaço-tempo curvo de Schwarzschild.  ║\n";
    std::cout << "╚══════════════════════════════════════════════════════════════╝\n";
    std::cout << "\n";
    
    std::cout << "Pressione ENTER para usar valores padrão ou digite novo valor:\n\n";
    
    std::string entrada;
    
    // Massa
    std::cout << "Massa do buraco negro [" << config.massa_solar << " M☉]: ";
    std::getline(std::cin, entrada);
    if (!entrada.empty()) {
        config.massa_solar = std::stod(entrada);
    }
    
    // Inclinação
    std::cout << "Ângulo de inclinação [" << config.angulo_inclinacao << "°]: ";
    std::getline(std::cin, entrada);
    if (!entrada.empty()) {
        config.angulo_inclinacao = std::stod(entrada);
    }
    
    // Resolução
    std::cout << "Largura da imagem [" << config.largura << " px]: ";
    std::getline(std::cin, entrada);
    if (!entrada.empty()) {
        config.largura = std::stoi(entrada);
    }
    
    std::cout << "Altura da imagem [" << config.altura << " px]: ";
    std::getline(std::cin, entrada);
    if (!entrada.empty()) {
        config.altura = std::stoi(entrada);
    }
    
    // Threads
    std::cout << "Número de threads [" << config.num_threads << "]: ";
    std::getline(std::cin, entrada);
    if (!entrada.empty()) {
        config.num_threads = std::stoi(entrada);
    }
}

// ============================================================
// FUNÇÃO PRINCIPAL
// ============================================================

int main(int argc, char* argv[]) {
    ConfiguracaoSimulacao config;
    bool apenas_analise = false;
    bool modo_interativo = false;
    
    // Opções de linha de comando
    static struct option opcoes_longas[] = {
        {"massa",       required_argument, nullptr, 'm'},
        {"spin",        required_argument, nullptr, 's'},
        {"inclinacao",  required_argument, nullptr, 'i'},
        {"distancia",   required_argument, nullptr, 'd'},
        {"largura",     required_argument, nullptr, 'W'},
        {"altura",      required_argument, nullptr, 'H'},
        {"fov",         required_argument, nullptr, 'f'},
        {"threads",     required_argument, nullptr, 't'},
        {"saida",       required_argument, nullptr, 'o'},
        {"analise",     no_argument,       nullptr, 'a'},
        {"interativo",  no_argument,       nullptr, 'I'},
        {"ajuda",       no_argument,       nullptr, 'h'},
        {nullptr, 0, nullptr, 0}
    };
    
    int opt;
    while ((opt = getopt_long(argc, argv, "m:s:i:d:W:H:f:t:o:aIh", opcoes_longas, nullptr)) != -1) {
        switch (opt) {
            case 'm':
                config.massa_solar = std::stod(optarg);
                break;
            case 's':
                config.spin = std::stod(optarg);
                if (config.spin > 0.01) {
                    config.tipo = TipoBuracoNegro::KERR;
                }
                break;
            case 'i':
                config.angulo_inclinacao = std::stod(optarg);
                break;
            case 'd':
                config.distancia_observador = std::stod(optarg);
                break;
            case 'W':
                config.largura = std::stoi(optarg);
                break;
            case 'H':
                config.altura = std::stoi(optarg);
                break;
            case 'f':
                config.fov = std::stod(optarg);
                break;
            case 't':
                config.num_threads = std::stoi(optarg);
                break;
            case 'o':
                config.diretorio_saida = optarg;
                break;
            case 'a':
                apenas_analise = true;
                break;
            case 'I':
                modo_interativo = true;
                break;
            case 'h':
                imprimir_ajuda(argv[0]);
                return 0;
            default:
                imprimir_ajuda(argv[0]);
                return 1;
        }
    }
    
    // Modo interativo se nenhum argumento
    if (argc == 1) {
        modo_interativo = true;
    }
    
    if (modo_interativo) {
        menu_interativo(config);
    }
    
    // Cria simulador
    SimuladorBuracoNegro simulador(config);
    
    // Apenas análise
    if (apenas_analise) {
        simulador.analise_fisica();
        return 0;
    }
    
    // Renderização
    if (!simulador.renderizar()) {
        std::cerr << "Erro na renderização!" << std::endl;
        return 1;
    }
    
    // Opção para análise adicional
    std::cout << "\nDeseja ver a análise física detalhada? (s/n): ";
    std::string resp;
    std::getline(std::cin, resp);
    if (resp == "s" || resp == "S" || resp == "sim") {
        simulador.analise_fisica();
    }
    
    std::cout << "\n✓ Simulação concluída com sucesso!\n\n";
    
    return 0;
}
