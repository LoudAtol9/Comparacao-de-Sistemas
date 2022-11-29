import numpy
import copy
import math

class Sistema :

    def __init__(self, matriz_sistema, matriz_resposta, bits=64) :

        self.matriz = copy.deepcopy(matriz_sistema)        # matriz com as constantes das variaveis
        self.matriz_resp = copy.deepcopy(matriz_resposta)  # matriz com as igualdades dos sistemas        
        self.matriz_exp = self.make_matriz_exp()           # as duas matrizes acima juntas

        self.TAM = len(self.matriz)                              # altura da matriz, largura = altura + 1
        self.resultado_k = [0] * (self.TAM)                      # lista com as variaveis
        self.met_pivotagem_parcial()                             # resolve o sistema com o maximo de precisao
        self.matriz_exp = self.make_matriz_exp()                 # reinicia o valor da matriz
        self.resultado_preciso = copy.deepcopy(self.resultado_k) # lista com o resultado preciso
        self.resultado_k = [0] * (self.TAM)                      # reinicia o valor da matriz
        self.orientacao_k = list(range(0, self.TAM))             # lembra o sentido da lista
        self.erro = [0] * (self.TAM)                             # erro real

        self.bits = bits     # num de bits usados  
        self.it_max = 1000   # num max de iteracoes 
        self.it_atual = 0    # contador de iteracoes
        self.E_tol = 0.0001 # precisao, tolerancia

        self.norm_nxn()

    def padrao_de_fabrica(self) :
        self.matriz_exp = self.make_matriz_exp()   # as duas matrizes acima juntas

        self.TAM = len(self.matriz_exp) - 1                 # altura da matriz, largura = altura + 1
        self.resultado_k = [numpy.float16(0)] * (self.TAM)  # lista com as variaveis
        self.orientacao_k = list(range(0, self.TAM))        # lembra o sentido da lista
        self.erro = 0                                       # erro real

        self.it_max = 1000   # num max de iteracoes 
        self.it_atual = 0    # contador de iteracoes
        self.E_tol = 0.00001 # precisao, tolerancia
        

    def calculo_erro_relativo(self) :
        for i in range(self.TAM) :
                if (self.resultado_k[i] != 0):
                    self.erro[i] = abs(self.resultado_preciso[i] - self.resultado_k[i])/abs(self.resultado_k[i])
                else:
                    self.erro[i] = abs(self.resultado_preciso[i] - self.resultado_k[i])


    # Recebe e devolve o número normalizado
    def norm_n(self, n):

        if self.bits == 16 :
            return numpy.float16(n)
        elif self.bits == 32 :
            return numpy.float32(n)
        elif self.bits == 64 :
            return numpy.float64(n)

    # Recebe e devolve uma matriz normalizada
    def norm_nxn(self):
        if self.bits == 16 :
            return self.normalizacao_16bits()
        elif self.bits == 32 :
            return self.normalizacao_32bits()
        elif self.bits == 64 :
            return self.normalizacao_64bits()

    # Todos os elementos da matriz em float16
    def normalizacao_16bits(self) :
        for i in range(len(self.matriz_exp)) :
            for j in range(len(self.matriz_exp) + 1) :
                self.matriz_exp[i][j] = numpy.float16(self.matriz_exp[i][j])

    # Todos os elementos da matriz em float32
    def normalizacao_32bits(self) :
        for i in range(len(self.matriz_exp)) :
            for j in range(len(self.matriz_exp) + 1) :
                self.matriz_exp[i][j] = numpy.float32(self.matriz_exp[i][j])
    
    # Todos os elementos da matriz em float64
    def normalizacao_64bits(self) :
        for i in range(len(self.matriz_exp)) :
            for j in range(len(self.matriz_exp) + 1) :
                self.matriz_exp[i][j] = numpy.float64(self.matriz_exp[i][j])
    
    # Retorna um número em notacao cientifica
    def notacao_cientifica(self, n) :
        return "{:4e}".format(n)
    

    # Cria a Matriz Expandida
    def make_matriz_exp(self) :

        self.matriz_exp = copy.deepcopy(self.matriz)

        for i in range(0, len(self.matriz)) :
            self.matriz_exp[i].append(self.matriz_resp[i])

        return self.matriz_exp



    # Metodo da Eliminacao de Gauss 
    def met_gauss_eliminacao(self) :

        try:

            # Resolve o sistema
            for j in range((self.TAM + 1) - 2) :

                # cria e armazena o pivo na casa correspondente a iteracao
                for i in range(j + 1,(self.TAM)):
                    self.matriz_exp[i][j] = (self.matriz_exp[i][j])/(self.matriz_exp[j][j])

                    # multiplica cada item da linha acima pelo pivo e
                    # a subtrai da linha de baixo
                    for z in range(j + 1,(self.TAM + 1)) :
                        self.matriz_exp[i][z] = self.matriz_exp[i][z] - self.matriz_exp[j][z] * self.matriz_exp[i][j]

                # passa pra próxima coluna e desce uma linha        
                j += 1

            # Acha o valor das constantes
            for u in range(self.TAM - 1, -1, -1) :

                # substitui o valor das constantes já achadas
                for v in range(self.TAM - 1, u, -1) :
                    self.matriz_exp[u][self.TAM] -= (self.matriz_exp[u][v] * self.resultado_k[v])

                # resolve a equacao
                self.resultado_k[u] = self.matriz_exp[u][self.TAM]/self.matriz_exp[u][u]

        except:
            self.resultado_k = [1] * self.TAM



    # Ordena em ordem decrescente a matriz pelas linhas
    def matriz_sort_linha(self, inicio) :

        for i in range(inicio, self.TAM):
            max = i
    
            for j in range(i + 1, self.TAM):
                # seleciona a linha maxima pra cada iteracao
                for z in range(inicio, self.TAM + 1) :
                    if abs(self.matriz_exp[j][z]) < abs(self.matriz_exp[max][z]):
                        break
                    if abs(self.matriz_exp[j][z]) > abs(self.matriz_exp[max][z]):
                        max = j
                        break

            # troca as linhas das matrizes
            (self.matriz_exp[i], self.matriz_exp[max]) = (copy.deepcopy(self.matriz_exp[max]), copy.deepcopy(self.matriz_exp[i]))




    # Ordena em ordem decrescente a matriz pelas colunas
    def matriz_sort_coluna(self, inicio) :

        for i in range(inicio, self.TAM):
            max = i
    
            for j in range(i + 1, self.TAM):
                # seleciona a coluna maxima pra cada iteracao
                for z in range(inicio, self.TAM) :
                    if abs(self.matriz_exp[z][j]) < abs(self.matriz_exp[z][max]):
                        break
                    if abs(self.matriz_exp[z][j]) > abs(self.matriz_exp[z][max]):
                        max = j
                        break

            # troca as colunas das matrizes
            for z in range(self.TAM):
                self.matriz_exp[z][i], self.matriz_exp[z][max] = self.matriz_exp[z][max], self.matriz_exp[z][i]
                self.orientacao_k[i] , self.orientacao_k[max] = self.orientacao_k[max] , self.orientacao_k[i]


    # Pega o maior elemento
    def matriz_sort_completa(self, inicio) :

        max_x = inicio
        max_y = inicio 

        for j in range(inicio + 1, self.TAM):

            # seleciona a coluna maxima pra cada iteracao
            for z in range(inicio, self.TAM) :
                if abs(self.matriz_exp[z][j]) > abs(self.matriz_exp[max_y][max_x]):
                    max_x = j
                    max_y = z
                        

        # troca as colunas das matrizes
        for z in range(self.TAM):
            self.matriz_exp[z][inicio], self.matriz_exp[z][max_x] = self.matriz_exp[z][max_x], self.matriz_exp[z][inicio]
        self.orientacao_k[inicio] , self.orientacao_k[max_x] = self.orientacao_k[max_x] , self.orientacao_k[inicio]
                

    # Function to find the partition position
    def partition(self, array, linha, low, high):
    
        # choose the rightmost element as pivot
        pivot = array[high][linha]
    
        # pointer for greater element
        i = low - 1
    
        # traverse through all elements
        # compare each element with pivot
        try:
            for j in range(low, high):
                if array[j][linha] >= self.modulo(pivot):
                
                    # If element smaller than pivot is found
                    # swap it with the greater element pointed by i
                    i = i + 1

                    # Swapping element at i with element at j
                    (array[i], array[j]) = (copy.deepcopy(array[j]), copy.deepcopy(array[i]))
        except TypeError:
            return float("NaN")
        # Swap the pivot element with the greater element specified by i
        (array[i + 1], array[high]) = (copy.deepcopy(array[high]), copy.deepcopy(array[i + 1]))
    
        # Return the position from where partition is done
        return i + 1
    
    # function to perform quicksort 
    def quickSort_2linhas(self, linha, low, high):
        if low < high:
        
            # Find pivot element such that
            # element smaller than pivot are on the left
            # element greater than pivot are on the right
            pi = self.partition(self.matriz_exp, linha, low, high)
            if math.isnan(pi): return float
    
            # Recursive call on the left of pivot
            self.quickSort_2linhas(linha, low, pi - 1)
    
            # Recursive call on the right of pivot
            self.quickSort_2linhas(linha, pi + 1, high)



    # Eliminacao de Gauss com Pivotagem Parcial
    def met_pivotagem_parcial(self) :

        # Resolve o sistema
        for j in range((self.TAM + 1) - 2) :

            # Pivotagem Parcial
            self.matriz_sort_linha(j)

            # cria e armazena o pivo na casa correspondente a iteracao
            for i in range(j + 1,(self.TAM)):

                self.matriz_exp[i][j] = (self.matriz_exp[i][j])/(self.matriz_exp[j][j])

                # multiplica cada item da linha acima pelo pivo e
                # a subtrai da linha de baixo
                for z in range(j + 1,(self.TAM + 1)) :
                    self.matriz_exp[i][z] = self.matriz_exp[i][z] - self.matriz_exp[j][z] * self.matriz_exp[i][j]

            # passa pra próxima coluna e desce uma linha        
            j += 1

        # Acha o valor das constantes
        for u in range(self.TAM - 1, -1, -1) :

            # substitui o valor das constantes já achadas
            for v in range(self.TAM - 1, u, -1) :
                self.matriz_exp[u][self.TAM] -= (self.matriz_exp[u][v] * self.resultado_k[v])

            # resolve a equacao
            self.resultado_k[u] = self.matriz_exp[u][self.TAM]/self.matriz_exp[u][u]



    # Eliminacao de Gauss com Pivotagem Completa
    def met_pivotagem_completa(self) :

        # Resolve o sistema
        for j in range((self.TAM + 1) - 2) :

            # Pivotagem Parcial
            self.matriz_sort_completa(j)
            self.matriz_sort_linha(j)

            # cria e armazena o pivo na casa correspondente a iteracao
            for i in range(j + 1,(self.TAM)):

                self.matriz_exp[i][j] = (self.matriz_exp[i][j])/(self.matriz_exp[j][j])

                # multiplica cada item da linha acima pelo pivo e
                # a subtrai da linha de baixo
                for z in range(j + 1,(self.TAM + 1)) :
                    self.matriz_exp[i][z] = self.matriz_exp[i][z] - self.matriz_exp[j][z] * self.matriz_exp[i][j]

            # passa pra próxima coluna e desce uma linha        
            j += 1

       # Acha o valor das constantes
        for u in range(self.TAM - 1, -1, -1) :

            # substitui o valor das constantes já achadas
            for v in range(self.TAM - 1, u, -1) :
                self.matriz_exp[u][self.TAM] -= (self.matriz_exp[u][v] * self.resultado_k[v])

            # resolve a equacao 
            self.resultado_k[u] = self.matriz_exp[u][self.TAM]/self.matriz_exp[u][u]

        # reordena o resultado
        resultado_copy = copy.deepcopy(self.resultado_k)
        for i in range(self.TAM):
            self.resultado_k[self.orientacao_k[i]] = resultado_copy[i]
        

    # Retorna o módulo de um número
    def modulo(self, i):
        if i < 0:
            return -i
        if i >= 0:
            return i


    # Metodo Gauss-Seidel
    def met_gauss_seidel(self) :

        resultado_k1 = [self.norm_n(0)] * (self.TAM) # x(k+1)
        soma = self.norm_n(0) # soma de x(k+1) e x(k)
        suma = 0
        continuar = True

        # Ordena a Matriz Diagonal Dominante
        for i in range(self.TAM) :
            self.matriz_sort_linha(i)

        if not continuar:
            self.resultado_k = [float("NaN")] * self.TAM
            return

        # Itera enquanto estiver fora da margem de erro e menor que a qnt max
        while (continuar and self.it_atual < self.it_max):

            # Percorre a coluna e reinicia o valor da soma
            for i in range(self.TAM):
                soma = self.norm_n(0)

                # Percorre linha, se for resposta adiciona, se for normal subtrai
                for j in range(self.TAM + 1) :
                    if (j == self.TAM):
                        soma += self.matriz_exp[i][j]
                    elif(j<i): 
                        soma -= self.matriz_exp[i][j] * resultado_k1[j]
                    elif(j>i): 
                        soma -=  self.matriz_exp[i][j] * self.resultado_k[j]
                    
                resultado_k1[i] = self.norm_n(soma/self.matriz_exp[i][i])

            # Verifica se o valor esta dentro da margem de erro
            for i in range(self.TAM) :
                if math.isnan(resultado_k1[i]) or math.isinf(resultado_k1[i]):
                    self.resultado_k = [float("NaN")] * self.TAM
                    return 0
                if self.modulo(self.resultado_k[i] - resultado_k1[i]) < self.norm_n(self.E_tol):
                    continuar = False

            # Passa o valor de k1 para k, iniciando outra iteracao
            self.resultado_k = copy.deepcopy(resultado_k1)

            self.it_atual += 1

        # reordena o resultado caso tenha ordenado a matriz
        resultado_copy = copy.deepcopy(self.resultado_k)
        for i in range(self.TAM):
            self.resultado_k[i] = resultado_copy[self.orientacao_k[i]]
    

if __name__ == '__main__':

    matriz_teste = [[2, -1, 10, -1], [0, 3, -1, 8], [10, -1, 2, 0], [-1, 11, -1, 3]]
    matriz_resposta = [-11, 15, 6, 25]

    sist = Sistema(matriz_teste, matriz_resposta, 16)
    print(sist.matriz_exp)
    sist.met_gauss_seidel()
    print(sist.resultado_k)
    print(sist.resultado_preciso)