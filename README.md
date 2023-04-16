# HPC-PROJECT

### Directory

```
.
│
├── LICENSE
├── Makefile
├── scripts   <--- script per eseguire misurazioni e test
│   ├── combine_result.sh
│   ├── gen-test-cases.sh
│   ├── measure.sh
│   ├── notify.sh
│   ├── test.sh
│   ├── varianza.sh
│   └── weak_scaling.sh
└── src
    ├── cuda-sph.cu  <--- File cuda che consegno in definitiva
    ├── hpc.h
    ├── omp-sph.c    <--- File omp che consegno in definitiva
    ├── sph.c
    ├── sph-cuda.cu  <--- parallelizzazione banale (sezione 4.1)
    ├── sph-cuda-shared.cu   <--- uso della shared (sezione 4.3)
    ├── sph-cuda-SoA.cu      <--- passaggio AoS a SoA (sezione 4.2)
    ├── sph-omp-dynamic.c    <--- partizionamento dinamico (sezione 3.2)
    ├── sph-omp-one-thread-pool.c   <--- partizionamento statico 1 thread pool (sezione 3.1.1)
    ├── sph-omp-simd-SoA.c   <--- omp + simd (sezione 3.3)
    └── sph-simd.c   <--- versione simd (sezione 2.1)

```

### Compiling
Per compilare tutti i sorgenti.
```bash
make all
```

Per compilare solo i sorgenti per la consegna
```bash
make consegna
```

Tutti gli eseguibili verranno messi nella cartella `bin`
Per cancellare gli eseguibili
```bash
make clean
```
