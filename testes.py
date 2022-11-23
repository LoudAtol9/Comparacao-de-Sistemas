import pytest
import numpy
import math
import prettytable
import time
import random
import copy
import resolucao_sistemas   


class Tabela:

    def __init__(self):

        self.table = ''
    
    def lst_metodos(self, n):

        if n == 1:
            return "| eliminacao | "
        if n == 2:
            return "|  parcial   | "
        if n == 3:
            return "|  completa  | "
        if n == 4:
            return "|   jacob    | "


    def add_linhas_vazias(self, qnt):
        self.table += "\n" * qnt

    def add_separacao_principal(self):
        self.table += "+----------------+---------+--------------+------------+--------------+------------------------------------------------------------------+\n"

    def add_sepracao_secundaria(self):
        self.table += "+------------+--------------+---------------+---------------------+---------------+\n"


    def add_cabecalho_principal(self):

        self.add_separacao_principal()
        self.table += "|     Matriz     | Tamanho | Complexidade |   Metodo   |    Tempo     |                          Erro Relativo                           |\n"
        self.add_separacao_principal()


    def add_linha_principal(self, cont, tam, complexidade, tipo, timer, erro):
        self.table += "|   " + ("{:9}".format(cont)) + "    |    " + str(tam) + "    |      " + str(complexidade) + "       " + self.lst_metodos(tipo) + str("{:7e}".format(timer)) + " | " + str(erro) + " |\n"


    def add_cabecalho_secundario(self):

        self.add_sepracao_secundaria()
        self.table += "|   Metodo   | Media Tempo  |  Media Erro   | Ocorrencias de NaN  | Erro sem NaN  |\n"
        self.add_sepracao_secundaria()


    def add_linha_secundaria(self, tipo, tempo, erro, NaN, erro2):
        self.table += self.lst_metodos(tipo) + str("{:12e}".format(tempo)) + " | " + str("{:12e}".format(erro)) + "% |       " + str("{:5}".format(NaN)) + "         | " + str("{:12e}".format(erro2)) + "  |\n"


class TestClass:
    
    def __init__(self, iteracoes):

        self.tempo_media_elim = 0
        self.tempo_media_parc = 0
        self.tempo_media_comp = 0
        self.tempo_media_jacob = 0

        self.erro_max_media_elim = 0
        self.erro_max_media_parc = 0
        self.erro_max_media_comp = 0
        self.erro_max_media_jacob = 0

        self.erro_max_media_NaN_elim = 0
        self.erro_max_media_NaN_parc = 0
        self.erro_max_media_NaN_comp = 0
        self.erro_max_media_NaN_jacob = 0

        self.NaN_elim = 0
        self.NaN_parc = 0
        self.NaN_comp = 0
        self.NaN_jacob = 0

        self.iter = iteracoes



    def random_matrizes(self, complexidade, tam, diagonal_dominante):

        m_matriz = []         # matriz final
        m_linha = [0] * tam  # matriz por linha
        d_matriz = []         # determina a diagonal dominante
        d_linha = [0] * tam # linha da diagonal dominante
        soma = 0

        inresolvivel = True

        # Determina o range dos numeros sorteados 
        if complexidade == 0:
            m_range = [-1, 1]
        if complexidade == 1:
            m_range = [-100, 100]
        if complexidade == 2:
            m_range = [-500, 500]
        if complexidade == 3:
            m_range = [-1000, 1000]

        while inresolvivel:
            # Criar matriz pra fazer a diagonal dominante
            for i in range(tam):
                if diagonal_dominante:
                    for j in range(tam):
                        if i == j:
                            d_linha[i] = 1
                d_matriz.append(copy.deepcopy(d_linha))

            # Embaralhar a Matriz
            if diagonal_dominante:
                for i in range(tam-1):
                    index_1 = random.randint(0, tam - 1)
                    index_2 = random.randint(0, tam - 1)
                    d_matriz[index_1], d_matriz[index_2] = copy.deepcopy(d_matriz[index_2]), copy.deepcopy(d_matriz[index_1])


            # Criar e sortear cada item da Matriz Principal
            for j in range(tam):
                soma = 0
                                                            
                for i in range(tam):
                    if d_matriz[j][i] == 0:
                        #while(i == j and a[i] == 0):
                        m_linha[i] = random.random() * random.randint(m_range[0], m_range[1])
                        soma += m_linha[i]
                    else:
                        m_linha[i] = '1'

                if(diagonal_dominante):  # elemento da diagonal principal maior que a soma dos outros
                    for i in range(tam):
                        if m_linha[i] == '1':
                            m_linha[i] = soma + (random.random() * random.randint(m_range[0], m_range[1]))/2

                m_matriz.append(copy.deepcopy(m_linha))

            # Sortear Matriz Resposta
            for i in range(tam):
                m_linha[i] = random.random() * random.randint(m_range[0], m_range[1])

            if(numpy.linalg.det(m_matriz) == 0):
                inresolvivel = True
            else:
                inresolvivel = False

        return m_matriz, m_linha





    def test_one(self, table, test_input, test_input_2, complex, cont=0, bits=16):

        sist = resolucao_sistemas.Sistema(test_input, test_input_2, bits)

        timer = time.perf_counter()
        sist.met_gauss_eliminacao()
        timer = time.perf_counter() - timer

        sist.calculo_erro_relativo()

        self.tempo_media_elim += timer/self.iter
        # Medida pra não dar erro 
        if math.isnan(max(sist.erro)/self.iter):
            self.NaN_elim += 1
            self.erro_max_media_elim += 1/self.iter
        else:
            self.erro_max_media_elim += max(sist.erro)/self.iter
            self.erro_max_media_NaN_elim += max(sist.erro)

        for i in range(sist.TAM):
            sist.erro[i] = sist.notacao_cientifica(sist.erro[i])

        table.add_linha_principal(cont, sist.TAM, complex, 1, timer, sist.erro)



    def test_two(self, table, test_input, test_input_2, complex, cont=0, bits=16):

        sist = resolucao_sistemas.Sistema(test_input, test_input_2, bits)

        timer = time.perf_counter()
        sist.met_pivotagem_parcial()
        timer = time.perf_counter() - timer

        sist.calculo_erro_relativo()

        self.tempo_media_parc += timer/self.iter
        if math.isnan(max(sist.erro)/self.iter):
            self.NaN_parc += 1
            self.erro_max_media_parc += 1/self.iter
        else:
            self.erro_max_media_parc += max(sist.erro)/self.iter
            self.erro_max_media_NaN_parc += max(sist.erro)

        for i in range(sist.TAM):
            sist.erro[i] = sist.notacao_cientifica(sist.erro[i])

        table.add_linha_principal(cont, sist.TAM, complex, 2, timer, sist.erro)
        


    def test_three(self, table, test_input, test_input_2, complex, cont=0, bits=16):

        sist = resolucao_sistemas.Sistema(test_input, test_input_2, bits)

        timer = time.perf_counter()
        sist.met_pivotagem_completa()
        timer = time.perf_counter() - timer

        sist.calculo_erro_relativo()

        self.tempo_media_comp += timer/self.iter
        if math.isnan(max(sist.erro)/self.iter):
            self.NaN_comp += 1
            self.erro_max_media_comp += 1/self.iter
        else:
            self.erro_max_media_comp += max(sist.erro)/self.iter
            self.erro_max_media_NaN_comp += max(sist.erro)

        for i in range(sist.TAM):
            sist.erro[i] = sist.notacao_cientifica(sist.erro[i])

        table.add_linha_principal(cont, sist.TAM, complex, 3, timer, sist.erro)



    def test_four(self, table, test_input, test_input_2, complex, cont=0, bits=16):

        sist = resolucao_sistemas.Sistema(test_input, test_input_2, bits)

        timer = time.perf_counter()
        sist.met_gauss_jacob()
        timer = time.perf_counter() - timer

        sist.calculo_erro_relativo()

        self.tempo_media_jacob += timer/self.iter
        if math.isnan(max(sist.erro)/self.iter):
            self.NaN_jacob += 1
            self.erro_max_media_jacob += 1/self.iter
        else:
            self.erro_max_media_jacob += max(sist.erro)/self.iter
            self.erro_max_media_NaN_jacob += max(sist.erro)

        for i in range(sist.TAM):
            sist.erro[i] = sist.notacao_cientifica(sist.erro[i])

        table.add_linha_principal(cont, sist.TAM, complex, 4, timer, sist.erro)






if __name__ == '__main__':

    iteracoes = int(input("Digite o número de Iteracoes: "))
    complexidade = int(input("Digite o range dos números sorteados (0 - 3): "))
    tam = int(input("Digite o tamanho da Matriz: "))
    bits = int(input("Digite em qnts Bits será a precisão (16 - 32 - 64): "))

    y = TestClass(iteracoes)
    table = Tabela()

    table.add_cabecalho_principal()
    
    for i in range(iteracoes):
        a, a1 = y.random_matrizes(complexidade, tam, False)
        y.test_one(table, a, a1, complexidade, i + 1, bits)
        y.test_two(table, a, a1,complexidade, i + 1, bits)
        y.test_three(table, a, a1, complexidade, i + 1, bits)
        y.test_four(table, a, a1, complexidade, i + 1, bits)
    
    table.add_separacao_principal()
    table.add_linhas_vazias(3)

    table.add_cabecalho_secundario()

    table.add_linha_secundaria(1, y.tempo_media_elim, y.erro_max_media_elim, y.NaN_elim, y.erro_max_media_NaN_elim/(iteracoes - y.NaN_elim))
    table.add_linha_secundaria(2, y.tempo_media_parc, y.erro_max_media_parc, y.NaN_parc, y.erro_max_media_NaN_parc/(iteracoes - y.NaN_parc))
    table.add_linha_secundaria(3, y.tempo_media_comp, y.erro_max_media_comp, y.NaN_comp, y.erro_max_media_NaN_comp/(iteracoes - y.NaN_comp))
    table.add_linha_secundaria(4, y.tempo_media_jacob, y.erro_max_media_jacob, y.NaN_jacob, y.erro_max_media_NaN_jacob/(iteracoes - y.NaN_jacob))

    table.add_sepracao_secundaria()


    f = open("outup.txt", 'w')
    f.write(table.table)

    f.close

