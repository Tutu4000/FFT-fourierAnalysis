import cmath
import matplotlib.pyplot as plt
import numpy as np

#Espectros visíveis
# Vermelho 480 a 405
# Laranja 510 a 480
# Amarelo 530 a 510
# Verde 600 a 530
# Ciano 620 a 600
# Azul 680 a 620
# Violeta 790 a 680

verde = np.linspace(530, 600, 100*7)
amarelo = np.linspace(510, 530, 100*5)
vermelho = np.linspace(405, 480, 100*3)
verde = np.sin(verde)
amarelo = np.sin(amarelo)
vermelho = np.sin(vermelho)

sinal = np.concatenate((verde, amarelo, vermelho))
plt.plot(sinal)
plt.show()
print(sinal)
def dft(lista):
    n = len(lista)
    result = []
    Z = cmath.exp(-2 * cmath.pi * 1j / n) #J é a unidade imaginária
    for k in range(n):
        soma = 0
        for t in range(n):
            soma += lista[t] * Z ** (t * k)
        result.append(soma)
    return result


def calcular_freq_domain(sinal):
    dft_sinal = dft(sinal)
    magnitudes = np.abs(dft_sinal) #magnitude = valor absoluto do número complexo
    T = len(sinal)
    frequencias = [i * 1/T for i in range(T)] #frequencia para um dado ponto i é i/T onde T é o tamanho do sinal. Isso vem da definição de frequencia que é 1/T
    magnitudes = magnitudes[:T//2] #magnitude e frequencia são simétricas (pq complexos tem conjugados) então pegamos só a metade
    frequencias = frequencias[:T//2]
    return frequencias, magnitudes


frequencias, magnitudes = calcular_freq_domain(sinal)
plt.plot(frequencias, magnitudes) #x = frequencias, y = magnitudes
#É como se magnitude fosse a força da frequencia
plt.title("Frequency map")
plt.show()

freqverde, magnitudesverde = calcular_freq_domain(verde)
plt.plot(freqverde, magnitudesverde, 'g', label='verde')

freqvermelho, magnitudesvermelho = calcular_freq_domain(vermelho)
plt.plot(freqvermelho, magnitudesvermelho, 'r', label='vermelho')

freqamarelo, magnitudesamarelo = calcular_freq_domain(amarelo)
plt.plot(freqamarelo, magnitudesamarelo, 'y', label='amarelo')


plt.title("Frequency map: Cores separadas vs Sinal completo")
plt.show()
