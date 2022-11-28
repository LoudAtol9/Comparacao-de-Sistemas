# Comparacao-de-Sistemas
Trabalho de Cálculo Numérico

## *Descrição* <h2>
O objetivo do programa é simular 4 métodos númericos e comparar suas carecteristicas entre si.<br> 
Os métodos são:

### *Eliminação de Gauss* <h4>
* Complexidade O( n^3 )<br>

A eliminação de Gauss, ou método de escalonamento, é um algoritmo para se resolver sistemas de equações lineares. Este método consiste em aplicar sucessivas operações elementares num sistema linear, para o transformar num sistema de mais fácil resolução, que apresenta exatamente as mesmas soluções que o original. Alta chance de ocorrer divisão com zero e não gerar resultado.<br>
<https://en.wikipedia.org/wiki/Gaussian_elimination>

### *Pivotagem Parcial* <h4>
* Complexidade O( n^3 ) + O( n log( n ) ) usando o Merge Sort<br>

Usa o mesmo racício que a eliminação de Gauss, porém a antes de gerar um novo pivô as linhas são ordenadas parar obter o maior elemento e conseguentemente reduzir o erro de ponto flutuante. Mais preciso que a Eliminação de Gauss, porém consome mais processamento e mais tempo.<br>
<https://en.wikipedia.org/wiki/Pivot_element>

### *Pivotagem Completa* <h4>
* Complexidade O( n^3 ) + O( n log( n ) ) usando o Merge Sort <br>

Usa o mesmo racício que a eliminação de Gauss e a Pivotagem parcial, antes de gerar um novo pivô as linhas e colunas são ordenadas parar obter o maior elemento dividindo no pivo conseguentemente reduzir o erro de ponto flutuante. Mais preciso que a Pivotagem Parical, porém consome um pouco mais de processamento e de tempo.<br>
<https://en.wikipedia.org/wiki/Pivot_element#Partial,_rook,_and_complete_pivoting>

### *Gauss - Seidel* <h4>
* Complexidade O( n^2 ) para cada iteração, sendo n o número de elementos na Matriz que não são zero <br>

Usa um racício difernete que a eliminação de Gauss, é um método iterativo, ou seja dificilmente chegará no resultado exato, mas conseguirá chegar bem próximo em pouca quantidade de tempo. Não funciona em todas as matrizes, ela deve ser diagonal dominante para convergir em algum resultado.<br>
<https://en.wikipedia.org/wiki/Gauss%E2%80%93Seidel_method><br><br>




## *Tutorial* <h2>

#### 1. *Requisitos Básicos* <h4>

* Python 3.8
* Instalador PIP

#### 2. *Bibliotecas Usadas* <h4>

* numpy
* random
* time
* math
* copy

Caso não tenha uma das Bibliotecas basta abrir o terminal e executar:
- $pip install "nome da biblioteca"

#### 3. *Download* <h4>
Clone o repositório usando o comando "$git clone https://github.com/LoudAtol9/Comparacao-de-Sistemas.git" ou baixe os arquivos pelo próprio github, ambos arquivos tem que estar na mesma pasta para o programa funcionar

#### 4. *Execução* <h4>
Basta executar o arquivo "Teste.py" colocar os parâmetros de teste ( n° de iterações, range, tamanho da matriz, n° de bits ) e após aguardar pouco, ou muito, tempo na mesma pasta será gerado um arquivo "output.txt" com todos os dados tabelados e no final do arquivo haverá as médias dos dados. 


