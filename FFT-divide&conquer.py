import cmath
import numpy as np
import matplotlib.pyplot as plt
import time
from numba import jit


def FFT(xn):
    N = len(xn)
    if N == 1:
        return xn
    else:
        Z = np.exp(-2j * np.pi * np.arange(N) / N)
        par = FFT(xn[::2])
        impar = FFT(xn[1::2])
        Aeven = par + Z[:int(N / 2)] * impar
        Aodd = par + Z[int(N / 2):] * impar
        return np.concatenate([Aeven, Aodd])

def DFT(lista):
    n = len(lista)
    result = []
    Z = cmath.exp(-2 * cmath.pi * 1j / n)  # J e a unidade imaginaria
    for k in range(n):
        soma = 0
        for t in range(n):
            soma += lista[t] * Z ** (t * k)
        result.append(soma)
    return result


def calcular_freq_domainDFT(sinal):
    dft_sinal = DFT(sinal)
    magnitudes = np.abs(dft_sinal)  # magnitude = valor absoluto do numero complexo
    T = len(sinal)
    frequencias = [i * 1 / T for i in range(
        T)]  # frequencia para um dado ponto i e i/T onde T e o tamanho do sinal. Isso vem da definicao de frequencia que e 1/T
    magnitudes = magnitudes[
                 :T // 2]  # magnitude e frequencia sao simetricas (pq complexos tem conjugados) entao pegamos so a metade
    frequencias = frequencias[:T // 2]
    return frequencias, magnitudes


def calcular_freq_domainFFT(sinal):
    fft_sinal = FFT(sinal)
    magnitudes = np.abs(fft_sinal)  # magnitude = valor absoluto do numero complexo
    T = len(sinal)
    frequencias = [i * 1 / T for i in range(
        T)]  # frequencia para um dado ponto i e i/T onde T e o tamanho do sinal. Isso vem da definicao de frequencia que e 1/T
    magnitudes = magnitudes[
                 :T // 2]  # magnitude e frequencia sao simetricas (pq complexos tem conjugados) entao pegamos so a metade
    frequencias = frequencias[:T // 2]
    return frequencias, magnitudes


def teste_um(tamanho):
    print("Teste para uma chamada:")
    verde = np.linspace(530, 600, 2 ** tamanho)
    verde = np.sin(verde)
    start = time.time()
    frequencias, magnitudes = calcular_freq_domainDFT(verde)
    end = time.time()
    decorrido = end - start
    print("DFT: {}".format(decorrido))
    plt.plot(frequencias, magnitudes, 'g', label='verde')
    start = time.time()
    frequencias, magnitudes = calcular_freq_domainFFT(verde)
    end = time.time()
    decorrido = end - start
    print("FFT: {}".format(decorrido))
    plt.plot(frequencias, magnitudes, 'b', label='verde')
    plt.title("Frequency map")
    plt.show()
    print("\n\n")


def teste_comparar_crescimento(size):
    print("Teste para comparar crescimento:")
    # função que vai aumentar gradualmente o tamanho do sinal e comparar o tempo de execução dos dois algoritmos
    # para cada tamanho de sinal
    tempos_dft = []
    tempos_fft = []
    tamanhos = [2 ** i for i in range(1, size)]
    for tamanho in tamanhos:
        sinal = np.linspace(530, 600, tamanho)
        sinal = np.sin(sinal)
        start = time.time()
        calcular_freq_domainDFT(sinal)
        end = time.time()
        decorrido = end - start
        tempos_dft.append(decorrido)
        start = time.time()
        calcular_freq_domainFFT(sinal)
        end = time.time()
        decorrido = end - start
        tempos_fft.append(decorrido)
    print("Tempos DFT: " + str(tempos_dft))
    print("Tempos FFT: " + str(tempos_fft))
    plt.plot(tamanhos, tempos_dft, 'g', label='DFT')
    plt.plot(tamanhos, tempos_fft, 'b', label='FFT')
    plt.title("Comparação de crescimento")
    plt.show()
    print("\n\n")


teste_comparar_crescimento(16)
teste_um(5)  # Argumento n onde tamanho = 2^n
