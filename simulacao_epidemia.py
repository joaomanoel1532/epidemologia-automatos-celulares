import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from scipy.integrate import odeint
from scipy.stats import linregress
from numba import jit

# MODELO SIR 

def sir_equations(y, t, N, beta, gamma):
    S, I, R = y
    dS_dt = -beta * S * I / N
    dI_dt = beta * S * I / N - gamma * I
    dR_dt = gamma * I
    
    return dS_dt, dI_dt, dR_dt

def run_sir_model(N, I0, R0, beta, gamma, total_time, time_steps):
    
    S0 = N - I0 - R0
    y0 = S0, I0, R0
    t = np.linspace(0, total_time, time_steps)
    solution = odeint(sir_equations, y0, t, args=(N, beta, gamma))
    S, I, R = solution.T
    
    return t, S, I, R

# AUTÔMATO CELULAR 

# Estados das Células
SUSCEPTIBLE = 0
INFECTED = 1
RECOVERED = 2

@jit(nopython=True)
def run_ca_simulation(L, time_steps, I_sir, gamma_ca, k_infection):
    
    # Inicializa a grade: todos suscetíveis, exceto um infectado no centro
    grid = np.full(L, SUSCEPTIBLE, dtype=np.int8)
    grid[L // 2] = INFECTED
    
    grid_history = np.zeros((time_steps, L), dtype=np.int8)
    grid_history[0, :] = grid
    
    for t in range(1, time_steps):
        new_grid = grid.copy()
        prevalence = I_sir[t] # Prevalência de infectados do SIR no tempo t
        
        for i in range(L):
            # Regra para células suscetíveis
            if grid[i] == SUSCEPTIBLE:
                # Verifica vizinhos
                left_neighbor = grid[(i - 1 + L) % L]
                right_neighbor = grid[(i + 1) % L]
                
                infected_neighbors = 0
                if left_neighbor == INFECTED:
                    infected_neighbors += 1
                if right_neighbor == INFECTED:
                    infected_neighbors += 1
                
                if infected_neighbors > 0:
                    # A probabilidade de infecção depende dos vizinhos e da prevalência geral
                    prob_infection = k_infection * infected_neighbors * prevalence
                    if np.random.rand() < prob_infection:
                        new_grid[i] = INFECTED
                        
            # Regra para células infectadas
            elif grid[i] == INFECTED:
                if np.random.rand() < gamma_ca:
                    new_grid[i] = RECOVERED
        
        grid = new_grid
        grid_history[t, :] = grid
        
    return grid_history

# ANÁLISE DE DIMENSÃO

def calculate_box_dimension(grid_state):
    # Apenas células não suscetíveis são contadas
    active_cells_indices = np.where(grid_state > 0)[0]
    
    if len(active_cells_indices) < 2: # Precisa de pelo menos 2 pontos para definir uma dimensão
        return 0.0

    L = len(grid_state)
    
    # Gera tamanhos de caixa
    max_power = int(np.log2(L))
    box_sizes = [2**i for i in range(max_power) if 2**i < L]
    
    if not box_sizes:
        return 0.0
        
    box_counts = []
    
    for size in box_sizes:
        count = 0
        for i in range(0, L, size):
            box = grid_state[i:i+size]
            if np.any(box > 0): # Verifica se há alguma célula ativa na caixa
                count += 1
        box_counts.append(count)

    # Evita log(0) e garante que temos pelo menos 2 pontos para a regressão
    valid_indices = [i for i, count in enumerate(box_counts) if count > 0]
    if len(valid_indices) < 2:
        return 0.0

    box_sizes_valid = [box_sizes[i] for i in valid_indices]
    box_counts_valid = [box_counts[i] for i in valid_indices]

    # Prepara dados para regressão linear: log(N(b)) vs log(1/b)
    log_box_sizes_inv = np.log([1.0/s for s in box_sizes_valid])
    log_counts = np.log(box_counts_valid)
    
    # Realiza a regressão linear para encontrar o coeficiente angular (a dimensão)
    try:
        slope, intercept, r_value, p_value, std_err = linregress(log_box_sizes_inv, log_counts)
    except ValueError:
        return 0.0
    
    return slope if not np.isnan(slope) else 0.0

# FUNÇÕES DE PLOTAGEM E EXECUÇÃO PRINCIPAL 

def plot_sir_curve(t, S, I, R, N, output_dir):
    """Gera o gráfico da curva epidêmica SIR."""
    plt.figure(figsize=(10, 6))
    plt.plot(t, S / N, 'b-', label='Suscetíveis (S)')
    plt.plot(t, I / N, 'r-', label='Infectados (I)')
    plt.plot(t, R / N, 'g-', label='Recuperados (R)')
    plt.xlabel('Tempo (dias)')
    plt.ylabel('Proporção da População')
    plt.title('Curva Epidêmica do Modelo SIR')
    plt.legend()
    plt.grid(True)
    plt.savefig(os.path.join(output_dir, "sir_curve.png"))
    plt.close()
    print("Gráfico da curva SIR salvo.")

def plot_ca_evolution(history, output_dir):
    """Gera o gráfico da evolução do autômato celular."""
    cmap = mcolors.ListedColormap(['#FFFFFF', '#FF0000', '#808080']) # Branco, Vermelho, Cinza
    bounds = [-0.5, 0.5, 1.5, 2.5]
    norm = mcolors.BoundaryNorm(bounds, cmap.N)

    plt.figure(figsize=(12, 8))
    plt.imshow(history, cmap=cmap, norm=norm, aspect='auto', interpolation='nearest')
    plt.xlabel('Posição na Grade')
    plt.ylabel('Tempo (dias)')
    plt.title('Evolução Espacial da Epidemia (Autômato Celular)')
    
    cbar = plt.colorbar(ticks=[0, 1, 2], orientation='vertical')
    cbar.set_ticklabels(['Suscetível', 'Infectado', 'Recuperado'])
    
    plt.savefig(os.path.join(output_dir, "ca_evolution.png"))
    plt.close()
    print("Gráfico da evolução do AC salvo.")

def plot_dimension_evolution(t, dimensions, output_dir):
    plt.figure(figsize=(10, 6))
    plt.plot(t, dimensions, 'm-', label='Dimensão de Similaridade (Dk)')
    plt.xlabel('Tempo (dias)')
    plt.ylabel('Dimensão de Similaridade')
    plt.title('Evolução da Complexidade Espacial da Epidemia')
    plt.ylim(bottom=0, top=max(dimensions)*1.1 if dimensions else 1)
    plt.legend()
    plt.grid(True)
    plt.savefig(os.path.join(output_dir, "dimension_evolution.png"))
    plt.close()
    print("Gráfico da dimensão de similaridade salvo.")

def main():
    # Parâmetros da Simulação 
    N = 100000        # População total
    I0 = 10           # Infectados iniciais
    R0_param = 0      # Recuperados iniciais
    beta = 0.5        # Taxa de transmissão
    gamma = 0.1       # Taxa de recuperação
    total_time = 150  # Tempo total de simulação em dias
    time_steps = 150  # Número de passos (1 passo por dia)
    
    # Parâmetros do Autômato Celular
    L = 512           # Tamanho da grade 1D
    gamma_ca = gamma  # Probabilidade de recuperação (consistente com SIR)
    k_infection = 0.8 # Fator de ajuste para probabilidade de infecção no AC

    # Diretório de saída
    output_dir = "data_output"
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    print("Executando o modelo SIR...")
    t, S, I, R = run_sir_model(N, I0, R0_param, beta, gamma, total_time, time_steps)
    
    print("Executando o Autômato Celular...")
    I_normalized = I / N # O AC usa a proporção de infectados
    ca_history = run_ca_simulation(L, time_steps, I_normalized, gamma_ca, k_infection)
    
    print("Calculando a dimensão de similaridade ao longo do tempo...")
    dimensions = [calculate_box_dimension(ca_history[i, :]) for i in range(time_steps)]

    # Geração dos Gráficos
    print("Gerando visualizações...")
    plot_sir_curve(t, S, I, R, N, output_dir)
    plot_ca_evolution(ca_history, output_dir)
    plot_dimension_evolution(t, dimensions, output_dir)

    print(f"\nSimulação concluída! Verifique os resultados na pasta '{output_dir}'.")

if __name__ == "__main__":
    main()