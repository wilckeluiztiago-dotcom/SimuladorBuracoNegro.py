#!/usr/bin/env python3
"""
SimuladorBuracoNegro - Animação de Geodésicas
Autor: Luiz Tiago Wilcke

Animação de partículas caindo em buraco negro
usando integração de geodésicas.
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter
from matplotlib.patches import Circle
import matplotlib.colors as mcolors

# ============================================================
# CONSTANTES
# ============================================================

G = 6.67430e-11
c = 299792458.0
c2 = c * c
M_SOL = 1.98892e30

def raio_schwarzschild(massa_solar):
    return 2.0 * G * massa_solar * M_SOL / c2

# ============================================================
# INTEGRADOR DE GEODÉSICAS
# ============================================================

class IntegradorGeodesica:
    """Integrador de geodésicas para Schwarzschild."""
    
    def __init__(self, massa_solar=10.0):
        self.M = massa_solar * M_SOL
        self.rs = raio_schwarzschild(massa_solar)
        
    def christoffel_r_tt(self, r):
        """Símbolo de Christoffel Γ^r_tt."""
        if r <= self.rs:
            return 0.0
        return self.rs * (r - self.rs) / (2.0 * r**3)
    
    def christoffel_r_rr(self, r):
        """Símbolo de Christoffel Γ^r_rr."""
        if r <= self.rs:
            return 0.0
        return -self.rs / (2.0 * r * (r - self.rs))
    
    def christoffel_r_phi_phi(self, r):
        """Símbolo de Christoffel Γ^r_φφ (plano equatorial)."""
        return -(r - self.rs)
    
    def christoffel_phi_r_phi(self, r):
        """Símbolo de Christoffel Γ^φ_rφ."""
        return 1.0 / r
    
    def derivadas(self, estado):
        """Calcula derivadas do estado [t, r, phi, u_t, u_r, u_phi]."""
        t, r, phi, u_t, u_r, u_phi = estado
        
        if r <= self.rs * 1.001:
            return np.zeros(6)
        
        # Derivadas das coordenadas
        dt = u_t
        dr = u_r
        dphi = u_phi
        
        # Derivadas das velocidades (equação geodésica)
        du_t = -2.0 * (self.rs / (2.0 * r * (r - self.rs))) * u_t * u_r
        du_r = (-self.christoffel_r_tt(r) * u_t**2 
                - self.christoffel_r_rr(r) * u_r**2
                - self.christoffel_r_phi_phi(r) * u_phi**2)
        du_phi = -2.0 * self.christoffel_phi_r_phi(r) * u_r * u_phi
        
        return np.array([dt, dr, dphi, du_t, du_r, du_phi])
    
    def passo_rk4(self, estado, h):
        """Passo Runge-Kutta 4ª ordem."""
        k1 = self.derivadas(estado)
        k2 = self.derivadas(estado + 0.5 * h * k1)
        k3 = self.derivadas(estado + 0.5 * h * k2)
        k4 = self.derivadas(estado + h * k3)
        
        return estado + h * (k1 + 2*k2 + 2*k3 + k4) / 6.0
    
    def integrar(self, r0, phi0, v_r0, L, num_passos=5000, h=0.001):
        """Integra trajetória de partícula massiva."""
        # Condições iniciais
        f = 1.0 - self.rs / r0
        u_t = 1.0 / f  # Energia = 1
        u_r = v_r0
        u_phi = L / (r0**2)
        
        estado = np.array([0.0, r0, phi0, u_t, u_r, u_phi])
        
        trajetoria = [estado.copy()]
        
        for _ in range(num_passos):
            estado = self.passo_rk4(estado, h)
            
            # Verifica horizonte
            if estado[1] <= self.rs * 1.001:
                break
            
            # Verifica escape
            if estado[1] > r0 * 10:
                break
                
            trajetoria.append(estado.copy())
        
        return np.array(trajetoria)

# ============================================================
# ANIMAÇÃO
# ============================================================

class AnimacaoBuracoNegro:
    """Animação de partículas caindo no buraco negro."""
    
    def __init__(self, massa_solar=10.0, num_particulas=8):
        self.massa_solar = massa_solar
        self.rs = raio_schwarzschild(massa_solar)
        self.integrador = IntegradorGeodesica(massa_solar)
        self.num_particulas = num_particulas
        
        # Cores para partículas
        self.cores = list(mcolors.TABLEAU_COLORS.values())
        
        # Gerar trajetórias
        self.trajetorias = []
        self.gerar_trajetorias()
        
    def gerar_trajetorias(self):
        """Gera trajetórias para todas as partículas."""
        print("Calculando geodésicas...")
        
        for i in range(self.num_particulas):
            # Condições iniciais variadas
            r0 = self.rs * (5 + i * 2)  # Distâncias diferentes
            phi0 = 2 * np.pi * i / self.num_particulas
            v_r0 = -0.1 * (0.5 + 0.5 * np.random.random())  # Velocidade radial para dentro
            L = np.sqrt(self.rs * r0) * (0.8 + 0.4 * np.random.random())  # Momento angular
            
            traj = self.integrador.integrar(r0, phi0, v_r0, L, num_passos=3000)
            self.trajetorias.append(traj)
            
        print(f"  {len(self.trajetorias)} trajetórias calculadas.")
    
    def para_cartesianas(self, traj):
        """Converte trajetória para coordenadas cartesianas."""
        r = traj[:, 1]
        phi = traj[:, 2]
        x = r * np.cos(phi)
        y = r * np.sin(phi)
        return x, y
    
    def criar_animacao(self, fps=30, duracao=10):
        """Cria animação."""
        fig, ax = plt.subplots(figsize=(10, 10))
        fig.patch.set_facecolor('black')
        ax.set_facecolor('black')
        
        # Limites
        limite = 30 * self.rs
        ax.set_xlim([-limite, limite])
        ax.set_ylim([-limite, limite])
        ax.set_aspect('equal')
        
        # Desenha buraco negro
        horizonte = Circle((0, 0), self.rs, color='black', zorder=10)
        ax.add_patch(horizonte)
        
        # Esfera de fótons
        fotons = Circle((0, 0), 1.5 * self.rs, fill=False, 
                        color='yellow', linestyle='--', alpha=0.5, zorder=5)
        ax.add_patch(fotons)
        
        # ISCO
        isco = Circle((0, 0), 3 * self.rs, fill=False,
                      color='green', linestyle=':', alpha=0.5, zorder=5)
        ax.add_patch(isco)
        
        # Disco de acreção (visual)
        for r_mult in np.linspace(3.5, 15, 12):
            disco = Circle((0, 0), r_mult * self.rs, fill=False,
                          color='orange', alpha=0.1 + 0.03 * (15 - r_mult), 
                          linewidth=5, zorder=1)
            ax.add_patch(disco)
        
        # Labels
        ax.set_xlabel('X (m)', color='white', fontsize=12)
        ax.set_ylabel('Y (m)', color='white', fontsize=12)
        ax.tick_params(colors='white')
        
        titulo = ax.set_title('', color='white', fontsize=14, pad=20)
        
        # Linhas de trajetória (rastro)
        linhas = []
        pontos = []
        
        for i, traj in enumerate(self.trajetorias):
            cor = self.cores[i % len(self.cores)]
            linha, = ax.plot([], [], color=cor, alpha=0.4, linewidth=1)
            ponto, = ax.plot([], [], 'o', color=cor, markersize=8)
            linhas.append(linha)
            pontos.append(ponto)
        
        # Prepara dados
        trajs_xy = [self.para_cartesianas(t) for t in self.trajetorias]
        max_frames = max(len(t) for t in self.trajetorias)
        num_frames = min(max_frames, fps * duracao)
        
        def init():
            for linha, ponto in zip(linhas, pontos):
                linha.set_data([], [])
                ponto.set_data([], [])
            titulo.set_text('')
            return linhas + pontos + [titulo]
        
        def animate(frame):
            tempo = frame / fps
            titulo.set_text(f'Buraco Negro {self.massa_solar:.0f} M☉  |  t = {tempo:.2f} s')
            
            for i, (linha, ponto, (x, y)) in enumerate(zip(linhas, pontos, trajs_xy)):
                # Índice proporcional
                idx = min(frame, len(x) - 1)
                
                # Rastro (últimos N pontos)
                inicio = max(0, idx - 100)
                linha.set_data(x[inicio:idx+1], y[inicio:idx+1])
                
                # Posição atual
                if idx < len(x):
                    ponto.set_data([x[idx]], [y[idx]])
                else:
                    ponto.set_data([], [])
            
            return linhas + pontos + [titulo]
        
        # Cria animação
        print("Gerando animação...")
        anim = FuncAnimation(fig, animate, init_func=init,
                           frames=num_frames, interval=1000/fps, blit=True)
        
        # Remove bordas
        for spine in ax.spines.values():
            spine.set_visible(False)
        
        ax.grid(True, alpha=0.1, color='gray')
        
        return anim, fig
    
    def salvar_gif(self, caminho, fps=30, duracao=10):
        """Salva animação como GIF."""
        anim, fig = self.criar_animacao(fps, duracao)
        
        print(f"Salvando GIF em {caminho}...")
        writer = PillowWriter(fps=fps)
        anim.save(caminho, writer=writer, dpi=100)
        plt.close(fig)
        
        print(f"✓ Animação salva: {caminho}")
    
    def mostrar(self, fps=30, duracao=10):
        """Mostra animação interativamente."""
        anim, fig = self.criar_animacao(fps, duracao)
        plt.show()

# ============================================================
# FUNÇÃO PRINCIPAL
# ============================================================

def main():
    print("╔══════════════════════════════════════════════════════════════╗")
    print("║     ANIMAÇÃO DE GEODÉSICAS - Buraco Negro                    ║")
    print("║     Autor: Luiz Tiago Wilcke                                 ║")
    print("╚══════════════════════════════════════════════════════════════╝")
    print()
    
    # Criar animação
    anim = AnimacaoBuracoNegro(massa_solar=10.0, num_particulas=6)
    
    # Salvar como GIF
    anim.salvar_gif('../saida/geodesicas.gif', fps=30, duracao=8)
    
    # Mostrar interativamente
    print("\nAbrindo animação interativa...")
    anim.mostrar(fps=30, duracao=8)

if __name__ == '__main__':
    main()
