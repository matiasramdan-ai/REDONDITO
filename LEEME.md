# REDONDITO

La `fftp` de SPECTER (`src/`), que se va convirtiendo en `sphp`: la FFT 2D de FFTW en (x,y) se reemplaza por la transformada en armónicos esféricos de SHTns en (θ,φ), plano por plano. Las rutinas escalares están en `src/sphp_sca.F90` (antes `fftp.F90`) y los módulos que comparten, `fftplans`, `threads` y `mpivars`, en `src/sphp_mod.F90` (antes `fftp_mod.F90`); las vectoriales, de (v_r,v_θ,v_φ) a sus coeficientes (Q,S,T) de SHTns y al revés (`fftp2d_vec_to_qst` y `fftp2d_qst_to_vec`), en `src/sphp_vec.F90`, que usa el plan y los cambios de reparto de `sphp_sca.F90`; los nombres de los módulos y de las rutinas todavía son los de SPECTER. Corre en paralelo (2026-10-31): las dos transformadas 2D hacen la transposición entre tareas, y la prueba da lo mismo con 1 a 6 tareas (`mpirun -np N bin/prueba_sphp`, con N hasta nz = 6). Se compila con `./compilar.sh`, que arma la biblioteca una vez y cuatro pruebas, cada una con su ejecutable en `bin/` (se corren con `mpirun -np N bin/<prueba>`): `prueba_sphp` (la 2D: nueve campos conocidos y una ida y vuelta al azar, N hasta 6), `prueba_vec` (las vectoriales: ocho campos con Q, S, T conocidos y una ida y vuelta al azar, N hasta 6), `prueba_z` (la FFT en z, periódica y con FC-Gram) y `prueba_3d` (las 3D, 2D más z, periódica y con FC-Gram). Las tablas de FC-Gram las leen de `SPECTER/tables`.

## Qué hace cada transformada con su entrada

Hay que distinguir dos cosas. Una transformada **en el lugar** deja el resultado en el mismo arreglo que recibe. Una que **destruye la entrada** deja el resultado en otro arreglo, pero usa la entrada como espacio de trabajo: después de la llamada la entrada no sirve y no se puede volver a usar. Estado al 2026-10-30:

| Rutina | Con qué | Entrada | Resultado |
|---|---|---|---|
| `fftp2d_real_to_complex_xy` | SHTns (`spat_to_SH`) | `in(nlat,nphi,z)` se **destruye**: SHTns hace la FFT en φ sobre el mismo arreglo (`INTENT(INOUT)`) | en otro arreglo, `out(lm,z)` |
| `fftp2d_complex_to_real_xy` | SHTns (`SH_to_spat`) | `in(lm,z)` **se conserva** (`INTENT(IN)`; medido: no cambia) | en otro arreglo, `out(nlat,nphi,z)` |
| `fftp2d_vec_to_qst` | SHTns (`spat_to_SHqst`) | `vr`, `vt`, `vp` se **destruyen** (medido, como en la escalar) | en otros arreglos, `q`, `s`, `t` |
| `fftp2d_qst_to_vec` | SHTns (`SHqst_to_spat`) | `q`, `s`, `t` **se conservan** (medido) | en otros arreglos, `vr`, `vt`, `vp` |
| `fftp1d_real_to_complex_z` | FFTW, con la continuación de FC-Gram | `inout(lm,z)` | **en el lugar**: entra f(lm,z) y sale f(lm,kz) |
| `fftp1d_complex_to_real_z` | FFTW | `inout(lm,kz)` | **en el lugar**: entra f(lm,kz) y sale f(lm,z) |
| `fftp3d_real_to_complex` | 2D + z | `in` se **destruye** (por la 2D) | en otro arreglo, `out(lm,kz)` |
| `fftp3d_complex_to_real` | z + 2D | `in` se **destruye** (la FFT en z es en el lugar) | en otro arreglo, `out(nlat,nphi,z)` |

En FFTW, la transformada directa (r2c) conserva su entrada y la inversa (c2r) la destruye; en SHTns es al revés: el análisis la destruye y la síntesis la conserva. Con los armónicos ortonormales, ninguna de las dos transformadas 2D hay que normalizarla.

## La transposición entre tareas

En el espacio físico cada tarea tiene sus planos, f(θ,φ,ksta:kend); en el espectral, todos los planos de sus lm, f(lmsta:lmend,nz). Es el reparto de SPECTER, que en el espectral reparte kx en lugar de lm. Pasar de un reparto al otro es la transposición entre tareas: el intercambio con `MPI_ISEND` y `MPI_IRECV`, que hacen `fftp3d_redist_k_to_lm` (de planos a lm, en la directa) y `fftp3d_redist_lm_to_k` (de lm a planos, en la inversa), las dos con los tipos de `fftp3d_create_block`, en el que cada tarea le manda a cada otra los lm de esa otra en sus propios planos (tipo `itype1`) y recibe de ella sus lm en los planos de la otra (tipo `itype2`). Esa no se puede evitar, porque cambia qué datos tiene cada tarea.

SPECTER hace además una transposición local, después de recibir en la directa y antes de mandar en la inversa: el ciclo por bloques `out(k,j,i) = c1(i,j,k)`. No cambia qué datos tiene la tarea sino el orden de los índices en memoria. La FFT 2D deja los datos como (kx,ky,z), con kx contiguo, y MPI los entrega en ese mismo orden; SPECTER los quiere como (z,ky,kx), con z contiguo, porque lo que sigue es la FFT en z a lo largo del primer índice.

En sphp no hace falta, porque el orden de los índices es el mismo antes y después del intercambio. SHTns deja cada plano como una columna de nlm coeficientes, `carr1(lm,z)`, con lm contiguo, y la salida es `out(lm,z)`, también con lm contiguo. Cada bloque que llega (mis lm en los planos de otra tarea) ocupa en `out` las columnas de esos planos con el mismo orden que tenía en `carr1`, así que MPI lo escribe directamente en su lugar. En la inversa pasa lo mismo al revés: se manda directo desde `in(lm,z)` y se recibe en `carr1(lm,z)`. Se ahorra una pasada entera sobre los datos con acceso salteado en memoria, que es lo que cuesta la transposición local. A cambio, la FFT en z (la que va a ser r) va a lo largo del segundo índice, con los datos salteados de a `lmend-lmsta+1`; FFTW lo resuelve con su plan para muchas transformadas a la vez.

## Validación (2026-10-31)

Las tres pruebas dan lo mismo con cualquier cantidad de tareas (probado de 1 a 6 la 2D y de 1 a 4 las otras dos). La 2D recupera los nueve campos conocidos con error de redondeo, entre 1e-30 y 1e-13. La FFT en z periódica da los coeficientes exactos de c cos(2πqz/L) + d sin(2πpz/L) a 3e-16, y la 3D, los de sin θ cos φ cos(2πqz/L) + cos θ sin(2πpz/L) + 1 a 9e-16. Con FC-Gram en z (C = 25, d = 5), la derivada en z de g(z) = sin(2.3z+0.4) + z² en [0,1], sola y multiplicada por sin θ cos φ + cos θ, tiene error relativo 2.5e-6 con nz = 64, 4.9e-8 con 128 y 2.2e-9 con 256 (orden 4.5 a 5.7), igual en z y en 3D; sin la continuación (control) el error sube a entre 18 y 108.

Comparación con la `fftp` del SPECTER original en 128³ (2026-11-01), z no periódica con FC-Gram (C = 25, d = 5), con la derivada en z de g(z) multiplicada por cos x + sin 2y en SPECTER y por sin θ cos φ + cos θ en sphp (nlat = nphi = 128, ℓmax = 63, el mayor que admite nphi = 128): las dos dan el mismo error, 4.94e-8 relativo, porque lo fija la continuación en z, que es la misma; la ida y vuelta da 9e-16 en SPECTER y 3e-14 en sphp (el redondeo de la transformada de Legendre). Tiempos por transformada 3D, directa o inversa, con una tarea: 16 ms en SPECTER y 11 ms en sphp; con 8, 3.3 y 1.5 ms. sphp es más rápida porque con ℓmax = 63 guarda 2080 coeficientes por plano contra 8320 de la FFT (128 × 65): la transposición mueve cuatro veces menos datos y la FFT en z hace cuatro veces menos columnas. Las pruebas de la comparación no quedaron en `src/`.

Normalización (decidida el 2026-10-31): todo viene normalizado, a diferencia de SPECTER, donde divide quien llama. La 2D de SHTns ya lo está, y la FFT en z directa (`fftp1d_real_to_complex_z`, y con ella la 3D directa) divide por nz después de FFTW, como el análisis de SHTns: los coeficientes son los del campo (cos(2πqz/L) da 1/2 en ±q) y la ida y vuelta en z o en 3D es la identidad. La inversa no normaliza.

Optimización polar de SHTns: eps = 1e-14, no la 1e-10 de MagIC (2026-11-01). Con 1e-10 la síntesis de modos de m alto pierde hasta 1e-11 (medido con una referencia en cuádruple precisión, ℓmax de 63 a 1023); con 1e-14 queda en el redondeo y tarda lo mismo. Detalles en `DISENO.md`.

Transformadas vectoriales (2026-11-01): `prueba_vec` da lo mismo con 1 a 6 tareas. Los campos radiales (Q), esferoidales (S) y toroidales (T, la rotación rígida entre ellos), de ℓ = 1 y de grado máximo, se recuperan con error de 1e-15 a 6e-14 en el análisis y de 1e-15 a 1e-14 en la síntesis; el toroidal sin(θ)^15 cos 15φ da 8e-13 en la síntesis, por el redondeo de la referencia (cos 15φ en el plano 6). La ida y vuelta al azar da 3e-14 a 9e-14: S y T redondean algo más que Q, porque las transformadas vectoriales multiplican y dividen por ℓ(ℓ+1).

