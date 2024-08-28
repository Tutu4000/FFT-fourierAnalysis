# Analysis of Algorithms - FFT and DFT

**Authors:**
- Arthur Faria Silva
- Geraldo Rodrigues de Melo Neto
- Maria Luísa Gabriel Domingues

**Title:** FFT and DFT - Examining Iterative and Divide and Conquer Techniques

---

## Introduction to the Chosen Problem

In the context of analyzing function values based on their frequency rather than their behavior over time, tools such as Fourier Series and Fourier Transform are widely used.

### Fourier Series

Fourier Series are used so that, in an infinite summation of sines and cosines, they can assume the behavior of the function being analyzed. With this generated series, it is possible to convert the function's interpretation domains: from the Time Domain, where the function's value is analyzed based on a parameter (usually time), to the Frequency Domain, where the most frequently assumed values by the function are analyzed within a parameter interval used in the Time Domain.

**Time Domain: Frequency Domain**

This conversion between domains is done by the Fourier Transform. Here we will address the Discrete Fourier Transform and its Divide and Conquer approach.

---

## Discrete Fourier Transform (DFT)

The Discrete Fourier Transform seeks to move the Observation Domain of a function by analyzing the vector behavior of a value present in the function. This is done by the summation:

\[ f_n = \sum_{n=0}^{N-1} f(n) \cdot e^{-\frac{2\pi i}{N}kn} \]

Where the function value in the Time Domain, called \(f\), has its value \(f_n\) (with \(n\) in the interval \(0..N-1\), and \(N\) being the number of “sample” values of \(f\)) subjected to a product with the correlation coefficient, composed of an exponential with an imaginary power (convertible to a pair \(\cos(n) + i\sin(n)\)) and with correlation degree \(k\), thus generating a value from the Complex Numerical Set, which is interpretable as a vector in a Cartesian plane Reals X Imaginaries.

This vector, in turn, in modulus, represents a part of the value \(F_k\) (resulting from the summation of samples for a certain degree \(k\)), which belongs to the Frequency Domain. The larger the modulus of this vector, the larger the result of \(F_k\).

---

## Derivation of DFT: Fast Fourier Transform (FFT)

The Fast Fourier Transform was designed to accelerate the calculation of the Discrete Fourier Transform by splitting the summation into two: even part (summation only with even \(n\) values) and odd part (summation only with odd \(n\) values), which can in turn be internally subdivided (considering the \(N/2\) elements as a new sequence \(0..N/2\)).

The possibility of applying the divide-and-conquer method in implementing this transform in code can be perceived, thus comparing the performance of DFT with its faster version, FFT.

---

## Comparative Analysis of the Algorithms

**DFT Algorithm:**

```python
def DFT(lista):
    n = len(lista)
    resultado = []
    z = e**(-2 * pi * i / n)
    for k in range(n):
        soma = 0
        for t in range(n):
            soma += lista[t] * z ** (t * k)
        resultado.append(soma)
    return resultado
```

**DFT(lista) Analysis:**

- Lines 2, 3, 4, and 10 execute only once, so we assign them a constant value \(c_1\).
- Assigning line 6 a cost \(c_2\), line 8 a cost \(c_3\), and line 9 a cost \(c_4\), we get that the algorithm's execution time function is:

\[ f(n) = an^2 + bn + c \]

It is possible to affirm that the algorithm has quadratic behavior:

\[ f(n) = Θ(g(n)) \]

If there exist positive constants \(c_1\), \(c_2\), and \(n_0\) such that:

\[ c_1g(n) \leq f(n) \leq c_2g(n) \]

As we want to demonstrate that \(f(n)\) has quadratic behavior:

\[ g(n) = n^2 \quad c_1n^2 \leq an^2 + bn + c \leq c_2n^2 \]

The inequality holds for \(c_1 = a/4\), \(c_2 = 7a/4\), and \(n_0 = 2 \cdot \max (|b|/a, \sqrt{|c|}/a)\). Therefore:

\[ f(n) = Θ(n^2) \]

**FFT Algorithm:**

```python
def FFT(N, lista):
    if N == 1:
        return lista
    else:
        lista_par = []
        lista_impar = []
        for i in range(N):
            if i % 2 == 0:
                lista_par.append(lista[i])
            else:
                lista_impar.append(lista[i])
        par = FFT(N / 2, lista_par)
        impar = FFT(N / 2, lista_impar)
        for k in range(N // 2):
            Aeven = par[k] + exp(-i * 2 * pi * k / N) * impar[k]
            Aodd = par[k + N // 2] + exp(-i * 2 * pi * (k + N // 2) / N) * impar[k + N // 2]
        return Aeven + Aodd
```

**FFT(N, lista) Analysis:**

The FFT algorithm is recursive and uses the divide and conquer technique. Thus, we specify:

- **Divide:** Split the problem into 2 subproblems (even[] and odd[]), where each one is approximately size \(n/2\).
- **Conquer:** Calculate the discrete Fourier transform recursively.
- **Combine:** Apply the two subsequences in the formula for the odd and even part of the Fourier transform and merge them at the end.

Thus, the execution time is expressed by recurrence, considering \(T(n)\), the execution time for an input of size \(n\), we should build a recurrence for \(T(n)\) in the form:

\[ T(n) = O(1) \quad \text{if } n \leq c \]

\[ T(n) = a \cdot T(n / b) + D(n) + C(n) \quad \text{if } n > c \]

Where \(c\) is a constant, \(D(n)\) is the Divide step, and \(C(n)\) is the Combine step.

If the list passed to the algorithm is of size 1, we just return the list. Thus, our base case is when we have \(c = 1\), with constant processing time \(O(1)\). Thus, FFT divides the problem into 2 subproblems of size \(n/2\), so the values of \(a\) and \(b\) are 2. The execution time is \(2T(n/2)\). For the divide step, we have a conditional loop to separate the values of even and odd positions, iterating over all \(N\) elements of the list, i.e., its processing time is \(\Theta(n)\). The combine step is done by a conditional loop of size \(n/2\), i.e., \(C(n) = n/2\), therefore it has a time equal to \(\Theta(n)\). Thus, the FFT recurrence is represented by:

\[ T(1) = O(1) \]

\[ T(n) = 2T(n/2) + n + n/2 \quad \text{if } n > 1 \]

**Case 2 of the Master Theorem:** If \(F(n) = \Theta(n \log_b a)\), then \(T(n) = \Theta(n \log_b a \cdot \log_b n)\)

As \(a = 2\), \(b = 2\), then \(n \log_b a = n \log_2 2 = n\). \(F(n) = 3n/2\). Then \(F(n) = 3n/2 = \Theta(n \log_2 2)\). Thus, by Case 2 of the Master Theorem, \(T(n) = \Theta(n \log_b a \cdot \log_b n)\), i.e., \(T(n) = \Theta(n \log n)\).

---

## Experimental Results of the Algorithms:

| Input Size | DFT Time (seconds) | FFT Time (seconds) |
|------------|--------------------|--------------------|
| 2^0        | 0.0                | 0.0                |
| 2^1        | 0.0                | 0.0                |
| 2^2        | 0.0                | 0.0                |
| 2^3        | 0.002              | 0.0                |
| 2^4        | 0.005              | 0.001              |
| 2^5        | 0.018              | 0.002              |
| 2^6        | 0.069              | 0.002             

 |
| 2^7        | 0.251              | 0.004              |
| 2^8        | 0.828              | 0.009              |
| 2^9        | 3.290              | 0.018              |
| 2^10       | 13.211             | 0.035              |
| 2^11       | 52.510             | 0.074              |

**X-axis:** Input Size, **Y-axis:** Time in seconds.  
*Note: The graph differs slightly from the table because a numba @jit compilation optimization was applied here.*

---

## Conclusions:

The Discrete Fourier Transform (DFT) is a powerful tool for frequency analysis in discrete signals. However, the Fast Fourier Transform (FFT) is an even more efficient implementation, which takes advantage of summation properties to split the larger summation into smaller, simpler calculations, thus drastically reducing computational complexity.

The \(O(n^2)\) time complexity of DFT can be prohibitive for large datasets. The improvement to \(O(n \log n)\) provided by FFT is significant and evident when observing empirical data.

Although FFT is optimized for inputs of size \(2^n\), it can be adapted to process inputs of any size. Moreover, for sufficiently large inputs, FFT will always be faster than the iterative version.

---

## References:

1. FFT. [S. l.], January 19, 2020. Available at: [link](https://fiscomp.if.ufrgs.br/index.php/FFT). Accessed on: February 16, 2024.
2. Understanding the Discrete Fourier Transform and the FFT. [S. l.: s. n.], 2024. Available at: [link](https://www.youtube.com/watch?v=QmgJmh2I3Fw). Accessed on: February 16, 2024.
3. Wavelets: a mathematical microscope. [S. l.: s. n.], 2023. Available at: [link](https://www.youtube.com/watch?v=jnxqHcObNK4). Accessed on: February 16, 2024.
4. Cooley, James; Tukey, John. An Algorithm for the Machine Calculation of Complex Fourier Series. [S. l.]. Available at: [link](https://web.stanford.edu/class/cme324/classics/cooley-tukey.pdf). Accessed on: February 16, 2024.
5. GitHub Repository: [https://github.com/Tutu4000/fourierAnalysis](https://github.com/Tutu4000/fourierAnalysis).

---
