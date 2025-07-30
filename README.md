# Modelo Híbrido SIR e Autômato Celular para Simulação de Epidemias

## 📌 Descrição do Projeto
Este projeto implementa um modelo híbrido para simular a dinâmica de doenças infecciosas, combinando a abordagem temporal do modelo **SIR (Suscetíveis, Infectados, Recuperados)** com a representação espacial de um **Autômato Celular (AC)**. O objetivo é não apenas modelar a curva epidêmica, mas também analisar a evolução da complexidade espacial do surto utilizando o conceito de **dimensão de similaridade**.

## 📊 Sobre os Modelos Utilizados
-   **Modelo SIR:** Um modelo compartimental clássico da epidemiologia matemática que descreve o fluxo de indivíduos entre os estados de Suscetível, Infectado e Recuperado ao longo do tempo. Ele fornece a dinâmica temporal da epidemia.
-   **Autômato Celular (AC):** Um modelo espacial discreto que consiste em uma grade de células. Aqui, ele é usado para simular a propagação da doença em um espaço unidimensional, onde o estado de cada célula (indivíduo) é atualizado com base em regras locais e na prevalência geral da doença fornecida pelo modelo SIR.
-   **Dimensão de Similaridade:** Um conceito da teoria dos fractais usado para quantificar a complexidade de um padrão espacial. Calculamos essa dimensão para os padrões gerados pelo AC para medir como a estrutura espacial da epidemia evolui.

## 🛠 Tecnologias Utilizadas
-   Python 3
-   NumPy
-   SciPy
-   Matplotlib
-   Numba (para otimização de performance)

## 📁 Estrutura do Projeto
```
📂 epdemologia
│-- 📂 data_output
│   └── ca_evoluton.png
│   └── dimension_evolution.png
│   └── sir_curve.png
│-- João Manoel - Uma Estrutura Computacional Híbrida Integrando a Dinâmica Temporal SIR com Autômatos Celulares Espaciais para Modelagem Epidemiológica.pdf
│-- 📜 simulacao_epidemia.py  # Script único contendo toda a implementação
│-- 📜 requirements.txt       # Dependências do projeto
│-- 📜 README.md              # Documentação do projeto
```
## 🔧 Como Executar o Projeto
1.  Clone este repositório:
    ```bash
    git clone [https://github.com/seu_usuario/Projeto_SIR_AC.git](https://github.com/seu_usuario/Projeto_SIR_AC.git)
    ```
2.  Acesse o diretório do projeto:
    ```bash
    cd epdemologia
    ```
3.  Instale as dependências necessárias:
    ```bash
    pip install -r requirements.txt
    ```
4.  Execute o script principal:
    ```bash
    python simulacao_epidemia.py
    ```

## 📈 Visualização dos Dados
O código gera uma série de visualizações que são salvas no diretório `data_output`:
1.  **Curva Epidêmica SIR:** Um gráfico mostrando a evolução das proporções de Suscetíveis, Infectados e Recuperados ao longo do tempo.
2.  **Evolução do Autômato Celular:** Uma imagem que mostra o estado da grade do AC em diferentes momentos da epidemia, visualizando a propagação espacial da doença.
3.  **Evolução da Dimensão de Similaridade:** Um gráfico que plota a dimensão de similaridade calculada em função do tempo, mostrando como a complexidade espacial do surto muda ao longo da epidemia.

## 📑 Fontes Teóricas
-   Os conceitos do modelo SIR foram baseados em "Modeling Infectious Diseases in Humans and Animals" de Keeling & Rohani.
-   Os conceitos de autômatos celulares e dimensão de similaridade foram baseados em "Cellular Automata: A Discrete View of the World" de Joel L. Schiff.
