# ============================================================
# SimuladorBuracoNegro - Sistema de Compilação
# Autor: Luiz Tiago Wilcke
# ============================================================

# Compilador e flags
CXX = g++
CXXFLAGS = -std=c++20 -O3 -march=native -Wall -Wextra -pthread
LDFLAGS = -pthread

# Debug flags (use: make DEBUG=1)
ifdef DEBUG
    CXXFLAGS = -std=c++20 -g -O0 -Wall -Wextra -pthread -DDEBUG
endif

# Diretórios
SRC_DIR = src
INC_DIR = include
BUILD_DIR = build
BIN_DIR = .

# Arquivos
SOURCES = $(wildcard $(SRC_DIR)/*.cpp)
OBJECTS = $(SOURCES:$(SRC_DIR)/%.cpp=$(BUILD_DIR)/%.o)
HEADERS = $(wildcard $(INC_DIR)/*.hpp)
TARGET = $(BIN_DIR)/simulador

# ============================================================
# ALVOS PRINCIPAIS
# ============================================================

.PHONY: all clean run debug help

all: $(TARGET)
	@echo ""
	@echo "╔══════════════════════════════════════════════════════════╗"
	@echo "║  Compilação concluída com sucesso!                       ║"
	@echo "║  Execute: ./simulador                                    ║"
	@echo "║  Para ajuda: ./simulador --ajuda                         ║"
	@echo "╚══════════════════════════════════════════════════════════╝"
	@echo ""

$(TARGET): $(OBJECTS)
	@echo "[LINK] Criando executável..."
	$(CXX) $(OBJECTS) -o $@ $(LDFLAGS)

$(BUILD_DIR)/%.o: $(SRC_DIR)/%.cpp $(HEADERS) | $(BUILD_DIR)
	@echo "[CXX] Compilando $<..."
	$(CXX) $(CXXFLAGS) -I$(INC_DIR) -c $< -o $@

$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)

# ============================================================
# LIMPEZA
# ============================================================

clean:
	@echo "[CLEAN] Removendo arquivos de compilação..."
	rm -rf $(BUILD_DIR)
	rm -f $(TARGET)
	@echo "Limpeza concluída."

# ============================================================
# EXECUÇÃO
# ============================================================

run: $(TARGET)
	./$(TARGET)

# Renderização rápida (baixa resolução)
quick: $(TARGET)
	./$(TARGET) -W 400 -H 300 -t 4

# Renderização de alta qualidade
hq: $(TARGET)
	./$(TARGET) -W 1920 -H 1080 -t 8

# Análise física apenas
analise: $(TARGET)
	./$(TARGET) --analise

# ============================================================
# DEBUG
# ============================================================

debug: clean
	$(MAKE) DEBUG=1

valgrind: debug
	valgrind --leak-check=full ./$(TARGET) -W 100 -H 100

# ============================================================
# INSTALAÇÃO
# ============================================================

install: $(TARGET)
	@echo "[INSTALL] Instalando em /usr/local/bin..."
	sudo cp $(TARGET) /usr/local/bin/simulador-buraco-negro
	@echo "Instalação concluída. Execute: simulador-buraco-negro"

uninstall:
	sudo rm -f /usr/local/bin/simulador-buraco-negro

# ============================================================
# AJUDA
# ============================================================

help:
	@echo ""
	@echo "SimuladorBuracoNegro - Comandos Disponíveis:"
	@echo ""
	@echo "  make          - Compila o projeto"
	@echo "  make run      - Compila e executa (modo interativo)"
	@echo "  make quick    - Renderização rápida (400x300)"
	@echo "  make hq       - Renderização alta qualidade (1920x1080)"
	@echo "  make analise  - Apenas análise física"
	@echo "  make clean    - Remove arquivos de compilação"
	@echo "  make debug    - Compila com símbolos de debug"
	@echo "  make install  - Instala no sistema"
	@echo "  make help     - Mostra esta mensagem"
	@echo ""
