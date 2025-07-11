# test# Tema 3 ASC

**Autor:** Artiom Pujleacov  
**Grupa:** 334CB

---

## Interpretare analiza Cachegrind

### Legenda

- **I refs:** numarul de referinte la instructiuni
- **I1/D1 misses:** ratari in cache-ul de nivel 1
- **D refs:** numar total de referinte la date (read/write)
- **LL misses:** ratari in cache-ul last level
- **Branches/Mispredicts:** salturi conditionale/indirecte si cat de des predictia a fost gresita

| Metrica         | `neopt` | `opt_m` | `blas`      |
|-----------------|---------|---------|-------------|
| I refs          | 229.794 | 229.794 | 4.248.922   |
| D refs          | 65.916  | 65.916  | 1.668.226   |
| D1 miss rate    | 4.9%    | 4.9%    | 2.5%        |
| LL miss rate    | 1.3%    | 1.3%    | 0.3%        |
| Branches        | 45.110  | 45.110  | 611.437     |
| Mispredict rate | 12.2%   | 12.2%   | 4.6%        |

### Interpretare

- **BLAS** are un numar foarte mare de instructiuni si referinte, dar o rata mica de ratari in cache-ul de nivel 1 si last level, ceea ce indica o buna utilizare a cache-ului. De asemenea, are un numar mare de salturi conditionale, dar cu o rata mica de predictie gresita, ceea ce sugereaza ca predictia salturilor este eficienta. BLAS are ramificari predictibile sau putine.
- In schimb, `neopt` si `opt_m` au un numar mult mai mic de instructiuni si referinte, dar o rata de ratari in cache-ul de nivel 1 si last level mai mare, ceea ce indica o utilizare mai putin eficienta a cache-ului. De asemenea, au un numar mai mic de salturi conditionale, dar o rata mai mare de predictie gresita, ceea ce sugereaza o predictie mai putin eficienta a salturilor. Rata mare de branch mispredict indica cod cu multe if-uri imprevizibile.

**Concluzie:**  
`blas` este cea mai eficienta implementare in ceea ce priveste utilizarea cache-ului si predictia salturilor, in timp ce `neopt` si `opt_m` au o utilizare mai putin eficienta a acestora.

---

## Efectul optimizarilor

Desi variantele `neopt` si `opt_m` au indicatori Cachegrind aproape identici (acelasi numar de referinte la instructiuni si date, aceleasi rate de ratari in cache), `opt_m` este semnificativ mai rapid in executia reala. Acest lucru se datoreaza optimizarilor facute manual in cod, care nu reduc neaparat numarul de accesari la memorie, dar:

- reduc numarul de operatii inutile (de exemplu prin eliminarea recalcularilor sau buclelor ineficiente)
- rearanjeaza accesarile la memorie pentru a imbunatati localitatea spatiala si temporala
- folosesc structuri de date mai eficiente (ex: bucle "interchange", blocuri, evitarea cache miss-urilor prin acces secvential etc.)

Cachegrind masoara doar numarul de accesari, dar nu surprinde complet beneficiile precum:

- rulare mai rapida pe CPU real datorita unui cod mai compact si mai bine vectorizat
- reducerea latentelor cauzate de accesari haotice
- eliminarea unor operatii complet prin optimizari (dead code, strength reduction, etc.)

Desi datele de cache sunt similare, `opt_m` are un timp de executie mai mic din cauza eficientizarii reale a codului si a utilizarii mai bune a procesorului.

---

## Analiza graficelor

### Timpii de executie

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

Am adaugat poza cu graficele in arhiva sub numele `timpi.png`.

### Observatii

- Implementarea `neopt` are cel mai mare timp de executie pentru valorile lui N. Timpul creste exponential cu dimensiunea, indicand un cod neoptimizat, cu acces ineficient la memorie si/sau fara vectorizare.
- Implementarea `opt_m` are un timp de executie mai mic decat `neopt`, dar mai mare decat `blas`. Timpul de executie creste exponential cu dimensiunea, dar nu la fel de rapid ca in cazul `neopt`, ceea ce sugereaza ca optimizarile au avut un impact semnificativ asupra performantei, dar nu la fel de mult ca optimizarile BLAS.
- Implementarea `blas` are cel mai mic timp de executie pentru toate dimensiunile, ceea ce sugereaza ca este cea mai eficienta implementare dintre cele trei. Timpul de executie creste mai lent decat in cazul `neopt` si `opt_m`, ceea ce sugereaza ca BLAS utilizeaza o abordare mai eficienta pentru a gestiona dimensiunile mai mari ale datelor. De asemenea, BLAS are un timp de executie constant pentru dimensiuni mai mari, ceea ce sugereaza ca este capabil sa gestioneze eficient dimensiuni mari ale datelor fara a afecta semnificativ performanta.

**Clasament eficienta:**  

1. BLAS  
2. opt_m  
3. neopt

---

## Explicatii implementari

### `neopt`

Aceasta versiune realizeaza calculele cerute folosind bucle clasice `for`, fara nicio forma de optimizare sau utilizare de biblioteci externe. Toate operatiile se fac in ordine:

1. Se calculeaza transpusa matricii A manual
2. Se face inmultirea matricii A transpuse cu vectorul x
3. Se inmulteste rezultatul cu matricea B
4. La final se adauga vectorul y

Acest cod e simplu si usor de inteles, dar nu este optimizat pentru performanta.

### `opt_m`

Aceasta versiune realizeaza aceleasi calcule, dar foloseste optimizari manuale:

- Reordonarea buclelor pentru a profita de cache-ul procesorului
- Evitarea calculelor repetate
- Utilizarea de blocuri sau prelucrari intermediare pentru a reduce dimensiunea datelor procesate la un moment dat

Aceasta versiune ofera o performanta buna fara a necesita biblioteci externe.

### `blas`

Aceasta versiune utilizeaza biblioteca BLAS (Basic Linear Algebra Subprograms), folosind `cblas_dgemv` pentru a efectua inmultirea matricii cu vectorul si `cblas_dgemm` pentru inmultirea matricii cu matrice.

---

## Prompt-uri si explicatii suplimentare

Am folosit ChatGPT 4.

Am folosit urmatoare succesiune de operatii pentru ca asa a fost in cerinta temei. Mai intai am facut blas-ul
si am modificat versiunea primita de la chat pana la una finala care nu foloseste calcule manuale. Apoi, am
inteles ca era mai usor sa o fac pe cea neoptimizata prima si dupa ea blas-ul si dupa optimizata. Dupa
ce am terminat cu blas-ul am facut versiunea neoptimizata, care a iesit foarte usor, pentru ca e simpla
ca pseudocod. Dupa aia am intrebat cum sa optimizez codul, si cu indiciile date de chat am ajuns la functia
finala
---

## Bucla N

Ultimele trei operatii din enunt implica lucrul cu vectori si matrici de dimensiune N, iar rezultatul final al fiecareia dintre aceste operatii este un vector coloana de dimensiune N x 1. Pentru a obtine acest vector, este necesar sa se calculeze fiecare element al sau in parte, deoarece fiecare pozitie a rezultatului depinde de valorile corespunzatoare din liniile matricilor implicate si din vectorii utilizati.

De aceea, se foloseste o bucla care se repeta de N ori, o data pentru fiecare linie a matricii. In cadrul fiecarei iteratii se realizeaza calculele specifice pentru pozitia respectiva din vectorul rezultat. Aceasta abordare permite construirea pas cu pas a vectorului final si este esentiala pentru a asigura corectitudinea si completitudinea rezultatului.

**In concluzie**, bucla de dimensiune N este esentiala pentru a obtine vectorul rezultat corect si complet, deoarece permite parcurgerea tuturor liniilor matricilor si vectorilor implicati in operatiile de inmultire si adunare.
