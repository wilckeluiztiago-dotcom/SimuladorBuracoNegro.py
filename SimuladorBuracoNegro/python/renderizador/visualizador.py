#!/usr/bin/env python3
"""
SimuladorBuracoNegro - Visualizador 3D
Autor: Luiz Tiago Wilcke

Visualização 3D de geodésicas e horizontes de buracos negros
usando Matplotlib.
"""

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.patches import Circle
import mpl_toolkits.mplot3d.art3d as art3d

# ============================================================
# CONSTANTES FÍSICAS
# ============================================================

G = 6.67430e-11          # Constante gravitacional (m³/(kg·s²))
c = 299792458.0          # Velocidade da luz (m/s)
c2 = c * c
M_SOL = 1.98892e30       # Massa do Sol (kg)

def raio_schwarzschild(massa_solar):
    """Calcula o raio de Schwarzschild para uma dada massa em massas solares."""
    M = massa_solar * M_SOL
    return 2.0 * G * M / c2

# ============================================================
# VISUALIZADOR 3D
# ============================================================

class VisualizadorBuracoNegro:
    """Visualizador 3D de buracos negros."""
    
    def __init__(self, massa_solar=10.0, spin=0.0):
        self.massa_solar = massa_solar
        self.spin = spin
        self.rs = raio_schwarzschild(massa_solar)
        
        # Cores
        self.cor_horizonte = 'black'
        self.cor_ergosfera = 'purple'
        self.cor_disco = cm.hot
        self.cor_geodesica = 'cyan'
        
    def criar_esfera(self, raio, resolucao=50):
        """Cria malha de esfera."""
        u = np.linspace(0, 2 * np.pi, resolucao)
        v = np.linspace(0, np.pi, resolucao)
        x = raio * np.outer(np.cos(u), np.sin(v))
        y = raio * np.outer(np.sin(u), np.sin(v))
        z = raio * np.outer(np.ones(np.size(u)), np.cos(v))
        return x, y, z
        
    def criar_disco(self, raio_interno, raio_externo, resolucao=100):
        """Cria malha de disco de acreção."""
        r = np.linspace(raio_interno, raio_externo, resolucao)
        theta = np.linspace(0, 2 * np.pi, resolucao)
        r, theta = np.meshgrid(r, theta)
        x = r * np.cos(theta)
        y = r * np.sin(theta)
        z = np.zeros_like(x)
        return x, y, z, r
    
    def temperatura_disco(self, r):
        """Perfil de temperatura Shakura-Sunyaev simplificado."""
        r_isco = 3.0 * self.rs
        if np.any(r < r_isco):
            r = np.maximum(r, r_isco)
        
        # T ∝ r^(-3/4) × [1 - (r_in/r)^(1/2)]^(1/4)
        x = r / r_isco
        T_max = 1e7  # Kelvin (exemplo)
        T = T_max * np.power(x, -0.75) * np.power(1.0 - np.sqrt(1.0 / x), 0.25)
        return T
    
    def desenhar_horizonte(self, ax, alpha=1.0):
        """Desenha o horizonte de eventos."""
        x, y, z = self.criar_esfera(self.rs)
        ax.plot_surface(x, y, z, color='black', alpha=alpha, shade=True)
        
    def desenhar_esfera_fotons(self, ax, alpha=0.3):
        """Desenha a esfera de fótons."""
        r_fotons = 1.5 * self.rs
        x, y, z = self.criar_esfera(r_fotons)
        ax.plot_surface(x, y, z, color='yellow', alpha=alpha, shade=False)
        
    def desenhar_disco_acrecao(self, ax, alpha=0.8):
        """Desenha o disco de acreção com mapa de temperatura."""
        r_isco = 3.0 * self.rs
        r_externo = 15.0 * self.rs
        
        x, y, z, r = self.criar_disco(r_isco, r_externo)
        T = self.temperatura_disco(r)
        T_norm = (T - T.min()) / (T.max() - T.min())
        
        ax.plot_surface(x, y, z, facecolors=self.cor_disco(T_norm), 
                       alpha=alpha, shade=False, linewidth=0)
    
    def desenhar_geodesica(self, ax, r0, phi0, theta0=np.pi/2, 
                           num_pontos=500, tipo='massiva'):
        """Desenha uma geodésica (trajetória de partícula)."""
        # Integração simplificada de geodésica
        dt = 0.1
        r = r0
        phi = phi0
        theta = theta0
        
        # Velocidade inicial (circular aproximada)
        v_phi = np.sqrt(G * self.massa_solar * M_SOL / r0) / r0
        v_r = 0.0
        
        trajetoria = []
        
        for _ in range(num_pontos):
            x = r * np.sin(theta) * np.cos(phi)
            y = r * np.sin(theta) * np.sin(phi)
            z = r * np.cos(theta)
            trajetoria.append([x, y, z])
            
            # Força gravitacional (Newtoniano para visualização)
            if r > self.rs * 1.01:
                a_r = -G * self.massa_solar * M_SOL / (r * r)
            else:
                break
            
            # Integração
            v_r += a_r * dt
            r += v_r * dt
            phi += v_phi * dt
            
            if r < self.rs * 1.01:
                break
        
        trajetoria = np.array(trajetoria)
        ax.plot3D(trajetoria[:, 0], trajetoria[:, 1], trajetoria[:, 2],
                  color=self.cor_geodesica, linewidth=1.5, alpha=0.8)
    
    def visualizar(self, mostrar_disco=True, mostrar_geodesicas=True,
                   mostrar_esfera_fotons=True, elevation=30, azimuth=45):
        """Cria visualização 3D completa."""
        fig = plt.figure(figsize=(14, 10))
        ax = fig.add_subplot(111, projection='3d')
        
        # Define limites
        limite = 20 * self.rs
        ax.set_xlim([-limite, limite])
        ax.set_ylim([-limite, limite])
        ax.set_zlim([-limite/2, limite/2])
        
        # Desenha componentes
        self.desenhar_horizonte(ax)
        
        if mostrar_esfera_fotons:
            self.desenhar_esfera_fotons(ax)
            
        if mostrar_disco:
            self.desenhar_disco_acrecao(ax)
            
        if mostrar_geodesicas:
            # Várias geodésicas em diferentes órbitas
            for r_mult in [5, 7, 10, 12]:
                for phi_inicio in [0, np.pi/2, np.pi, 3*np.pi/2]:
                    self.desenhar_geodesica(ax, r_mult * self.rs, phi_inicio)
        
        # Configurações visuais
        ax.set_facecolor('black')
        ax.xaxis.pane.fill = False
        ax.yaxis.pane.fill = False
        ax.zaxis.pane.fill = False
        ax.xaxis.pane.set_edgecolor('darkgray')
        ax.yaxis.pane.set_edgecolor('darkgray')
        ax.zaxis.pane.set_edgecolor('darkgray')
        
        ax.set_xlabel('X (m)', color='white')
        ax.set_ylabel('Y (m)', color='white')
        ax.set_zlabel('Z (m)', color='white')
        ax.tick_params(colors='white')
        
        ax.view_init(elev=elevation, azim=azimuth)
        
        # Título
        titulo = f'Buraco Negro de Schwarzschild\n'
        titulo += f'M = {self.massa_solar:.0f} M☉, rs = {self.rs/1000:.1f} km'
        ax.set_title(titulo, color='white', fontsize=14, pad=20)
        
        plt.tight_layout()
        return fig, ax
    
    def salvar(self, caminho, **kwargs):
        """Salva a visualização em arquivo."""
        fig, ax = self.visualizar(**kwargs)
        fig.savefig(caminho, dpi=150, facecolor='black', 
                   edgecolor='none', bbox_inches='tight')
        plt.close(fig)
        print(f"Visualização salva em: {caminho}")

# ============================================================
# VISUALIZADOR DE POTENCIAL EFETIVO
# ============================================================

def plotar_potencial_efetivo(massa_solar=10.0, L_valores=None):
    """Plota o potencial efetivo para diferentes momentos angulares."""
    rs = raio_schwarzschild(massa_solar)
    
    if L_valores is None:
        L_valores = [3.0, 3.5, 4.0, 4.5, 5.0]  # Em unidades de rs
    
    r = np.linspace(1.5 * rs, 30 * rs, 1000)
    
    fig, ax = plt.subplots(figsize=(12, 8))
    fig.patch.set_facecolor('#1a1a2e')
    ax.set_facecolor('#16213e')
    
    for L in L_valores:
        L_real = L * rs * c  # Momento angular real
        
        # Potencial efetivo: V²_eff = (1 - rs/r)(1 + L²/(r²c²))
        V_eff_sq = (1.0 - rs/r) * (1.0 + L_real**2 / (r**2 * c2))
        V_eff = np.sqrt(np.maximum(V_eff_sq, 0))
        
        ax.plot(r/rs, V_eff, label=f'L = {L:.1f} rs·c', linewidth=2)
    
    # Linhas de referência
    ax.axvline(x=1.0, color='red', linestyle='--', alpha=0.5, 
               label='Horizonte (r = rs)')
    ax.axvline(x=1.5, color='yellow', linestyle='--', alpha=0.5,
               label='Esfera de fótons (r = 1.5rs)')
    ax.axvline(x=3.0, color='green', linestyle='--', alpha=0.5,
               label='ISCO (r = 3rs)')
    
    ax.set_xlabel('r / rs', color='white', fontsize=12)
    ax.set_ylabel('Potencial Efetivo V_eff', color='white', fontsize=12)
    ax.set_title(f'Potencial Efetivo - Buraco Negro de {massa_solar:.0f} M☉', 
                 color='white', fontsize=14)
    
    ax.tick_params(colors='white')
    ax.legend(loc='upper right', facecolor='#16213e', 
              edgecolor='white', labelcolor='white')
    ax.grid(True, alpha=0.3, color='gray')
    
    ax.set_xlim([1, 30])
    ax.set_ylim([0.9, 1.1])
    
    for spine in ax.spines.values():
        spine.set_edgecolor('white')
    
    plt.tight_layout()
    return fig, ax

# ============================================================
# FUNÇÃO PRINCIPAL
# ============================================================

def main():
    """Função principal para demonstração."""
    print("╔══════════════════════════════════════════════════════════════╗")
    print("║     VISUALIZADOR DE BURACO NEGRO - Python                    ║")
    print("║     Autor: Luiz Tiago Wilcke                                 ║")
    print("╚══════════════════════════════════════════════════════════════╝")
    print()
    
    # Parâmetros
    massa = 10.0  # Massas solares
    
    # Criar visualizador
    viz = VisualizadorBuracoNegro(massa_solar=massa)
    
    # Informações
    print(f"Massa: {massa} M☉")
    print(f"Raio de Schwarzschild: {viz.rs:.2f} m ({viz.rs/1000:.2f} km)")
    print(f"ISCO: {3*viz.rs:.2f} m")
    print()
    
    # Visualização 3D
    print("[1/3] Gerando visualização 3D...")
    viz.salvar('../saida/visualizacao_3d.png', elevation=30, azimuth=45)
    
    # Potencial efetivo
    print("[2/3] Gerando gráfico de potencial efetivo...")
    fig, ax = plotar_potencial_efetivo(massa)
    fig.savefig('../saida/potencial_efetivo.png', dpi=150, 
                facecolor='#1a1a2e', bbox_inches='tight')
    plt.close(fig)
    print("Gráfico salvo em: ../saida/potencial_efetivo.png")
    
    # Visualização interativa
    print("[3/3] Abrindo visualização interativa...")
    viz.visualizar()
    plt.show()
    
    print("\n✓ Visualização concluída!")

if __name__ == '__main__':
    main()
