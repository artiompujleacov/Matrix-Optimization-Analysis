# test# Tema 3 ASC

**Autor:** Artiom Pujleacov  
**Grupa:** 334CB

---

## Interpretare analiză Cachegrind

### Legendă

- **I refs:** numărul de referințe la instrucțiuni
- **I1/D1 misses:** ratări în cache-ul de nivel 1
- **D refs:** număr total de referințe la date (read/write)
- **LL misses:** ratări în cache-ul last level
- **Branches/Mispredicts:** salturi condiționale/indirecte și cât de des predicția a fost greșită

| Metrica         | `neopt` | `opt_m` | `blas`      |
|-----------------|---------|---------|-------------|
| I refs          | 229.794 | 229.794 | 4.248.922   |
| D refs          | 65.916  | 65.916  | 1.668.226   |
| D1 miss rate    | 4.9%    | 4.9%    | 2.5%        |
| LL miss rate    | 1.3%    | 1.3%    | 0.3%        |
| Branches        | 45.110  | 45.110  | 611.437     |
| Mispredict rate | 12.2%   | 12.2%   | 4.6%        |

### Interpretare

- **BLAS** are un număr foarte mare de instrucțiuni și referințe, dar o rată mică de ratări în cache-ul de nivel 1 și last level, ceea ce indică o bună utilizare a cache-ului. De asemenea, are un număr mare de salturi condiționale, dar cu o rată mică de predicție greșită, ceea ce sugerează că predicția salturilor este eficientă. BLAS are ramificări predictibile sau puține.
- În schimb, `neopt` și `opt_m` au un număr mult mai mic de instrucțiuni și referințe, dar o rată de ratări în cache-ul de nivel 1 și last level mai mare, ceea ce indică o utilizare mai puțin eficientă a cache-ului. De asemenea, au un număr mai mic de salturi condiționale, dar o rată mai mare de predicție greșită, ceea ce sugerează o predicție mai puțin eficientă a salturilor. Rata mare de branch mispredict indică cod cu multe if-uri imprevizibile.

**Concluzie:**  
`blas` este cea mai eficientă implementare în ceea ce privește utilizarea cache-ului și predicția salturilor, în timp ce `neopt` și `opt_m` au o utilizare mai puțin eficientă a acestora.

---

## Efectul optimizărilor

Deși variantele `neopt` și `opt_m` au indicatori Cachegrind aproape identici (același număr de referințe la instrucțiuni și date, aceleași rate de ratări în cache), `opt_m` este semnificativ mai rapid în execuția reală. Acest lucru se datorează optimizărilor făcute manual în cod, care nu reduc neapărat numărul de accesări la memorie, dar:

- reduc numărul de operații inutile (de exemplu prin eliminarea recalculărilor sau buclelor ineficiente)
- rearanjează accesările la memorie pentru a îmbunătăți localitatea spațială și temporală
- folosesc structuri de date mai eficiente (ex: bucle "interchange", blocuri, evitarea cache miss-urilor prin acces secvențial etc.)

Cachegrind măsoară doar numărul de accesări, dar nu surprinde complet beneficiile precum:

- rulare mai rapidă pe CPU real datorită unui cod mai compact și mai bine vectorizat
- reducerea latențelor cauzate de accesări haotice
- eliminarea unor operații complet prin optimizări (dead code, strength reduction, etc.)

Deși datele de cache sunt similare, `opt_m` are un timp de execuție mai mic din cauza eficientizării reale a codului și a utilizării mai bune a procesorului.

---

## Analiza graficelor

### Timpii de execuție

#### `opt_m`

```
N=200:  Time=0.042094
N=400:  Time=0.322833
N=600:  Time=1.076220
N=800:  Time=2.536097
N=1000: Time=4.948999
```

#### `neopt`

```
N=200:  Time=0.241318
N=400:  Time=2.037641
N=600:  Time=7.376293
N=800:  Time=16.650749
N=1000: Time=28.030949
```

#### `blas`

```
N=200:  Time=0.008862
N=400:  Time=0.066758
N=600:  Time=0.190512
N=800:  Time=0.425241
N=1000: Time=0.839287
```

Am adăugat poza în arhivă sub numele `timpi.png`.

![Screenshot 2025-05-18 182627](https://github.com/user-attachments/assets/4ddb72a7-99b0-4143-afa6-1c6476c17da5)


### Observații

- Implementarea `neopt` are cel mai mare timp de execuție pentru valorile lui N. Timpul crește exponențial cu dimensiunea, indicând un cod neoptimizat, cu acces ineficient la memorie și/sau fără vectorizare.
- Implementarea `opt_m` are un timp de execuție mai mic decât `neopt`, dar mai mare decât `blas`. Timpul de execuție crește exponențial cu dimensiunea, dar nu la fel de rapid ca în cazul `neopt`, ceea ce sugerează că optimizările au avut un impact semnificativ asupra performanței, dar nu la fel de mult ca optimizările BLAS.
- Implementarea `blas` are cel mai mic timp de execuție pentru toate dimensiunile, ceea ce sugerează că este cea mai eficientă implementare dintre cele trei. Timpul de execuție crește mai lent decât în cazul `neopt` și `opt_m`, ceea ce sugerează că BLAS utilizează o abordare mai eficientă pentru a gestiona dimensiunile mai mari ale datelor. De asemenea, BLAS are un timp de execuție constant pentru dimensiuni mai mari, ceea ce sugerează că este capabil să gestioneze eficient dimensiuni mari ale datelor fără a afecta semnificativ performanța.

**Clasament eficiență:**  

1. BLAS  
2. opt_m  
3. neopt

---

## Explicații implementări

### `neopt`

Această versiune realizează calculele cerute folosind bucle clasice `for`, fără nicio formă de optimizare sau utilizare de biblioteci externe. Toate operațiile se fac în ordine:

1. Se calculează transpusa matricei A manual
2. Se face înmulțirea matricei A transpuse cu vectorul x
3. Se înmulțește rezultatul cu matricea B
4. La final se adaugă vectorul y

Acest cod e simplu și ușor de înțeles, dar nu este optimizat pentru performanță.

### `opt_m`

Această versiune realizează aceleași calcule, dar folosește optimizări manuale:

- Reordonarea buclelor pentru a profita de cache-ul procesorului
- Evitarea calculelor repetate
- Utilizarea de blocuri sau prelucrări intermediare pentru a reduce dimensiunea datelor procesate la un moment dat

Această versiune oferă o performanță bună fără a necesita biblioteci externe.

### `blas`

Această versiune utilizează biblioteca BLAS (Basic Linear Algebra Subprograms), folosind `cblas_dgemv` pentru a efectua înmulțirea matricei cu vectorul și `cblas_dgemm` pentru înmulțirea matricei cu matrice.

---

## Prompt-uri și explicații suplimentare

Am folosit ChatGPT 4 pentru a rafina implementările și a obține variantele optimizate și neoptimizate, precum și pentru a înțelege mai bine avantajele folosirii BLAS și a optimizărilor manuale.

```
request:
C = B * At
D = Ct * A
i = 0
while (i < N) {
    y = Ct * x
    x = C * y
    i += 1
}
y = D * x

hi,i need to implement these operations only using blas operations,this is the signature of the function:
double* my_solver(int N, double *A, double *B, double *x) {
	return NULL;
}



raspuns:
To implement the following operations using only BLAS routines inside the my_solver function:

c
C = B * At        // At is A transposed
D = Ct * A        // Ct is C transposed
loop N times:
    y = Ct * x
    x = C * y
y = D * x
You can use BLAS Level 3 routine cblas_dgemm for matrix-matrix multiplication and BLAS Level 2 routine cblas_dgemv for matrix-vector multiplication.

✅ Solution:
Here's how you can implement it:

c
#include <cblas.h>
#include <stdlib.h>

double* my_solver(int N, double *A, double *B, double *x) {
    // Allocate matrices and vectors
    double *C = (double *)malloc(N * N * sizeof(double));
    double *Ct = (double *)malloc(N * N * sizeof(double));
    double *D = (double *)malloc(N * N * sizeof(double));
    double *y = (double *)malloc(N * sizeof(double));
    double *tmp = (double *)malloc(N * sizeof(double));

    // C = B * At  -> At is A transposed
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans,
                N, N, N, 1.0, B, N, A, N, 0.0, C, N);

    // Ct = C^T  (we transpose C manually)
    for (int i = 0; i < N; ++i)
        for (int j = 0; j < N; ++j)
            Ct[j * N + i] = C[i * N + j];

    // D = Ct * A
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                N, N, N, 1.0, Ct, N, A, N, 0.0, D, N);

    // Copy x to tmp (x will be updated)
    for (int i = 0; i < N; ++i)
        tmp[i] = x[i];

    // Repeat N times:
    //   y = Ct * x
    //   x = C * y
    for (int iter = 0; iter < N; ++iter) {
        cblas_dgemv(CblasRowMajor, CblasNoTrans,
                    N, N, 1.0, Ct, N, tmp, 1, 0.0, y, 1);

        cblas_dgemv(CblasRowMajor, CblasNoTrans,
                    N, N, 1.0, C, N, y, 1, 0.0, tmp, 1);
    }

    // Final y = D * x (x is now in tmp)
    cblas_dgemv(CblasRowMajor, CblasNoTrans,
                N, N, 1.0, D, N, tmp, 1, 0.0, y, 1);

    // Cleanup
    free(C);
    free(Ct);
    free(D);
    free(tmp);

    return y;
}
```

```
request:
how can i implement the manual transpose of the matrix using blas functions

raspuns:
BLAS does not provide a built-in function to transpose a matrix manually. The BLAS interface is designed for high-performance linear algebra operations, and transposition is typically handled implicitly via arguments in BLAS routines (e.g., using CblasTrans in cblas_dgemm).

✅ Option 1: Use Transposition Flags in BLAS
Instead of manually transposing, you can tell BLAS to treat the matrix as transposed using the TransA or TransB flags.

Example:

c
cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, ...);  // Use Aᵗ directly
So you don't need to transpose A or C manually for:

c
C = B * At
→ use:

c
cblas_dgemm(..., B, ..., A, ..., CblasTrans);
And for:

c
D = Ct * A
→ use:

c
cblas_dgemm(..., C, ..., A, ..., CblasTrans, CblasNoTrans); // Treat C as transposed
✅ Option 2: Manual Transpose (only if required)
If you really need a separate transposed copy of a matrix, do it like this:

c
void transpose_matrix(int N, double *src, double *dst) {
    for (int i = 0; i < N; ++i)
        for (int j = 0; j < N; ++j)
            dst[j * N + i] = src[i * N + j]; // dst is the transpose of src
}
But again: this is outside of BLAS, and it's rarely necessary unless some other library requires the transposed matrix to be stored explicitly.
```

```
request:
now using the same signature,make a function just to manually do the operations i mentioned without using any external libraries

double* my_solver(int N, double *A, double *B, double *x) {
	return NULL;
}

raspuns:
Sure! Here's the complete function with a brief explanation and proper comments to help you understand what's happening at each step:

This implementation follows the operation sequence:

java
C = B * At       // At is the transpose of A
D = Ct * A       // Ct is the transpose of C
for i in 0 to N:
    y = Ct * x
    x = C * y
y = D * x        // Final result
Below is the full C function:

c
#include "utils.h"
#include <stdlib.h>
#include <string.h>

/*
 *
 *
 *
 *   C = B * At
 *   D = Ct * A
 *   for i = 0..N:
 *     y = Ct * x
 *     x = C * y
 *   y = D * x
 *
 */
double *my_solver(int N, double *A, double *B, double *x)
{
	double *At = (double *)malloc(N * N * sizeof(double));
	double *C = (double *)malloc(N * N * sizeof(double));
	double *Ct = (double *)malloc(N * N * sizeof(double));
	double *D = (double *)malloc(N * N * sizeof(double));
	double *y = (double *)malloc(N * sizeof(double));
	double *tmp = (double *)malloc(N * sizeof(double));
	double *result = (double *)malloc(N * sizeof(double));

	for (int i = 0; i < N; ++i)
		for (int j = 0; j < N; ++j)
			At[j * N + i] = A[i * N + j];

	// C = B * At
	for (int i = 0; i < N; ++i)
		for (int j = 0; j < N; ++j)
		{
			C[i * N + j] = 0.0;
			for (int k = 0; k < N; ++k)
				C[i * N + j] += B[i * N + k] * At[k * N + j];
		}

	// Ct = transpose C
	for (int i = 0; i < N; ++i)
		for (int j = 0; j < N; ++j)
			Ct[j * N + i] = C[i * N + j];

	// D = Ct * A
	for (int i = 0; i < N; ++i)
		for (int j = 0; j < N; ++j)
		{
			D[i * N + j] = 0.0;
			for (int k = 0; k < N; ++k)
				D[i * N + j] += Ct[i * N + k] * A[k * N + j];
		}

	memcpy(tmp, x, N * sizeof(double));

	//y = Ct * x; x = C * y;
	for (int iter = 0; iter < N; ++iter)
	{
		// y = Ct * tmp
		for (int i = 0; i < N; ++i)
		{
			y[i] = 0.0;
			for (int j = 0; j < N; ++j)
				y[i] += Ct[i * N + j] * tmp[j];
		}

		// tmp = C * y
		for (int i = 0; i < N; ++i)
		{
			tmp[i] = 0.0;
			for (int j = 0; j < N; ++j)
				tmp[i] += C[i * N + j] * y[j];
		}
	}

	// result = D * tmp
	for (int i = 0; i < N; ++i)
	{
		result[i] = 0.0;
		for (int j = 0; j < N; ++j)
			result[i] += D[i * N + j] * tmp[j];
	}

	free(At);
	free(C);
	free(Ct);
	free(D);
	free(y);
	free(tmp);

	return result;
}
✅ This implementation is pure C, without using BLAS or any external libraries.

💡 If performance becomes important, this same structure can later be rewritten using cblas_dgemm and cblas_dgemv from BLAS for significant speedups.

Let me know if you want that BLAS version too!
```

```
request:
now suggest a way to optimize the give function to make it faster,without using external 
librarier,without chnaging the complexity of the operations

raspuns:
To optimize your function without using external libraries and without changing the algorithmic complexity, you can apply cache-friendly and low-level optimizations. These focus on improving data locality, reducing memory access cost, and avoiding redundant computations.

✅ Suggested Optimizations:
1. Loop Order Optimization for Better Cache Locality
Matrix multiplications currently use:

c
for (i)
  for (j)
    for (k)
      result[i][j] += A[i][k] * B[k][j];
This accesses matrix B[k][j] column-wise (non-contiguous in memory), which is cache-inefficient. You can change the loop order to:

c
for (i)
  for (k)
    for (j)
      result[i][j] += A[i][k] * B[k][j];
This makes all three accesses mostly row-wise and improves CPU cache performance significantly.

2. Manual Blocking (Tiling)
Divide large matrices into small blocks to take better advantage of the CPU cache (L1/L2).

Example block size:

c
#define BLOCK_SIZE 32
Then modify the matrix multiplication to work in blocks:

c
for (ii = 0; ii < N; ii += BLOCK_SIZE)
  for (kk = 0; kk < N; kk += BLOCK_SIZE)
    for (jj = 0; jj < N; jj += BLOCK_SIZE)
      for (i = ii; i < min(ii + BLOCK_SIZE, N); ++i)
        for (k = kk; k < min(kk + BLOCK_SIZE, N); ++k)
          for (j = jj; j < min(jj + BLOCK_SIZE, N); ++j)
            C[i * N + j] += A[i * N + k] * B[k * N + j];
This improves spatial locality and cache reuse.

3. Avoid Recomputing Transposes
Instead of repeatedly computing the transpose, compute it once and reuse it. You already do this for At and Ct, which is great — just ensure you never re-transpose or access transposed matrices inefficiently afterward.

4. Use restrict Keyword (Optional)
If you are compiling with a modern compiler (e.g., GCC/Clang), you can tell the compiler that pointers don’t alias:

c
double *restrict A, *restrict B, *restrict C;
This allows the compiler to optimize more aggressively.

5. Minimize memcpy
The line:

c
memcpy(tmp, x, N * sizeof(double));
is fine, but if x is not reused afterward, consider just working with tmp = x directly (carefully, only if you're sure x won't be reused or overwritten unexpectedly). This saves memory copy time.

6. Unroll Small Loops
For small fixed-size loops like:

c
for (int j = 0; j < N; ++j)
If N is small and fixed, manually unroll the loop or let the compiler do it using #pragma unroll.

Summary Table
Optimization	Gains	Cost / Risk
Loop order change (ijk → ikj)	High (cache-friendly)	None
Manual tiling	Very High	More code
Avoid recomputing transposes	Moderate	Already done
Use restrict keyword	Moderate	Compiler-dependent
Minimize memcpy usage	Low	Code clarity
Manual loop unrolling	Low-Moderate	Code complexity
```

---

## Bucla N

Ultimele trei operații din enunț implică lucrul cu vectori și matrici de dimensiune N, iar rezultatul final al fiecăreia dintre aceste operații este un vector coloană de dimensiune N x 1. Pentru a obține acest vector, este necesar să se calculeze fiecare element al său în parte, deoarece fiecare poziție a rezultatului depinde de valorile corespunzătoare din liniile matricilor implicate și din vectorii utilizați.

De aceea, se folosește o buclă care se repetă de N ori, o dată pentru fiecare linie a matricii. În cadrul fiecărei iterații se realizează calculele specifice pentru poziția respectivă din vectorul rezultat. Această abordare permite construirea pas cu pas a vectorului final și este esențială pentru a asigura corectitudinea și completitudinea rezultatului.

**În concluzie**, bucla de dimensiune N este esențială pentru a obține vectorul rezultat corect și complet, deoarece permite parcurgerea tuturor liniilor matricilor și vectorilor implicați în operațiile de înmulțire și adunare.
