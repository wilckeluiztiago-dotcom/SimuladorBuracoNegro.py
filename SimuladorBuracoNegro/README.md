# Simulador de Buraco Negro Relativístico

![C++20](https://img.shields.io/badge/C++-20-blue.svg)
![Python 3](https://img.shields.io/badge/Python-3.8+-green.svg)
![License](https://img.shields.io/badge/License-MIT-yellow.svg)

**Autor:** Luiz Tiago Wilcke

Simulador de buracos negros com ray tracing relativístico em espaço-tempo curvo de Schwarzschild. Gera imagens fotorrealistas mostrando lente gravitacional, disco de acreção e sombra do buraco negro.

---

## 🌌 Características

- **Ray tracing relativístico** em espaço-tempo de Schwarzschild
- **Disco de acreção** com perfil de temperatura Shakura-Sunyaev
- **Efeitos físicos**: Doppler beaming, redshift gravitacional, lente gravitacional
- **Multithreading** para renderização rápida
- **Visualização Python** com Matplotlib (3D e animações)
- **Análise física** completa (temperatura de Hawking, entropia, etc.)

---

## 📁 Estrutura do Projeto

```
SimuladorBuracoNegro/
├── cpp/
│   ├── include/
│   │   ├── constantes.hpp      # Constantes físicas fundamentais
│   │   ├── schwarzschild.hpp   # Métrica de Schwarzschild
│   │   ├── kerr.hpp            # Métrica de Kerr (rotativo)
│   │   ├── integrador.hpp      # Integrador geodésico RK4
│   │   ├── disco_acrecao.hpp   # Modelo de disco de acreção
│   │   ├── ray_tracer.hpp      # Ray tracing relativístico
│   │   └── simulador.hpp       # Classe principal
│   ├── src/
│   │   └── main.cpp            # Programa principal
│   └── Makefile
├── python/
│   ├── renderizador/
│   │   └── visualizador.py     # Visualização 3D
│   └── animacao/
│       └── animacao_geodesica.py  # Animação de partículas
├── saida/                       # Imagens geradas
└── README.md
```

---

## 🚀 Compilação

### Requisitos
- GCC 10+ ou Clang 12+ (suporte C++20)
- Make
- Python 3.8+ com NumPy e Matplotlib (opcional)

### Compilar
```bash
cd cpp
make
```

### Executar
```bash
./simulador              # Modo interativo
./simulador --ajuda      # Ver opções
```

---

## 📖 Uso

### Linha de Comando

```bash
# Renderização padrão (800x600)
./simulador

# Personalizado
./simulador -m 20 -i 60 -W 1920 -H 1080 -t 8
#           │    │     │          │      └── 8 threads
#           │    │     └──────────┴── Full HD
#           │    └── 60° inclinação
#           └── 20 massas solares

# Apenas análise física
./simulador --analise -m 100
```

### Opções

| Opção | Descrição | Padrão |
|-------|-----------|--------|
| `-m, --massa` | Massa em M☉ | 10 |
| `-i, --inclinacao` | Ângulo de visão (graus) | 75 |
| `-W, --largura` | Largura da imagem (px) | 800 |
| `-H, --altura` | Altura da imagem (px) | 600 |
| `-t, --threads` | Número de threads | 4 |
| `-a, --analise` | Mostrar análise física | - |

### Atalhos Make

```bash
make quick    # 400x300 (teste rápido)
make hq       # 1920x1080 (alta qualidade)
make analise  # Apenas análise física
```

---

## 🔬 Física Implementada

### Métrica de Schwarzschild
```
ds² = -(1 - rs/r)dt² + (1 - rs/r)⁻¹dr² + r²dΩ²
```

### Raio de Schwarzschild
```
rs = 2GM/c² ≈ 2.95 km × (M/M☉)
```

### Temperatura de Hawking
```
T = ℏc³/(8πGMk) ≈ 6×10⁻⁸ K × (M☉/M)
```

### Disco de Acreção (Shakura-Sunyaev)
```
T(r) ∝ r⁻³/⁴ × [1 - (r_in/r)^½]^¼
```

---

## 🐍 Visualização Python

```bash
cd python/renderizador
python3 visualizador.py

cd ../animacao
python3 animacao_geodesica.py
```

---

## 📊 Performance

| Resolução | Tempo | Taxa |
|-----------|-------|------|
| 200×150 | 0.06s | 519K px/s |
| 800×600 | 0.30s | 1.6M px/s |
| 1920×1080 | ~1.3s | 1.6M px/s |

---

## 📄 Licença

MIT License - Luiz Tiago Wilcke © 2024

---

## 🔗 Referências

- Misner, Thorne & Wheeler - *Gravitation* (1973)
- Shakura & Sunyaev - *Disk Model* (1973)
- Event Horizon Telescope - *First Image* (2019)
