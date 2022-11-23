import numpy
import copy

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
        self.E_tol = 0.00001 # precisao, tolerancia

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




    # Metodo Gauss-Jacob
    def met_gauss_jacob(self) :

        resultado_k1 = [numpy.float16(0)] * (self.TAM) # x(k+1)
        soma = numpy.float16(0) # soma de x(k)

        continuar = True

        for i in range(self.TAM) :
            self.matriz_sort_linha(i)

        while (continuar and self.it_atual < self.it_max):
            for i in range(self.TAM) :
                soma = numpy.float16(0)

                for j in range(0, self.TAM) :
                    if j != i :
                        soma += self.matriz_exp[i][j] * self.resultado_k[j]

                resultado_k1[i] = self.norm_n(1/self.matriz_exp[i][i]) * (self.matriz_exp[i][self.TAM] - soma)


            self.resultado_k = copy.deepcopy(resultado_k1)

            self.it_atual += 1
        




    # Metodo Gauss-Seidel
    def met_gauss_seidel(self) :

        resultado_k1 = [0] * (self.TAM) # x(k+1)
        soma = 0 # soma de x(k+1) e x(k)

        continuar = True

        while (continuar and self.it_atual < self.it_max):
            for i in range(self.TAM):

                resultado_k1[i] = self.matriz_exp[self.TAM - 1][i]/self.matriz_exp[i][i]
                soma = 0

                for j in range(self.TAM) :
                    if(j<i): soma += self.matriz_exp[i][j] * resultado_k1[j]
                    if(j>i): soma +=  self.matriz_exp[i][j] * self.resultado_k[j]
                    
                resultado_k1[i] -= soma/self.matriz_exp[i][i]

            #for i in range(self.TAM) :
            #    if abs(self.resultado_k[i] - resultado_k1[i]) < self.E_tol:
            #        continuar = False

            self.resultado_k = copy.deepcopy(resultado_k1)

            self.it_atual += 1

        # reordena o resultado
        resultado_copy = copy.deepcopy(self.resultado_k)
        for i in range(self.TAM):
            self.resultado_k[i] = resultado_copy[self.orientacao_k[i]]
    
    
