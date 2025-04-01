"""
Módulo de alineamiento local de secuencias biológicas (Smith-Waterman).

Contiene funciones relacionadas con el cálculo de alineamientos entre pares de secuencias, así como la comparación de rendimiento con herramientas estándar como Biopython.

Funciones incluidas:
- `smith_waterman`: Implementación manual del algoritmo Smith-Waterman.
- `medir_rendimiento_custom`: Mide el tiempo de ejecución de la implementación propia.
- `medir_rendimiento_biopython`: Mide el tiempo usando PairwiseAligner de Biopython.
- `guardar_resultados_alineamientos`: Guarda alineamientos obtenidos.
- `cargar_matriz_sustitucion`: Carga matrices como BLOSUM o PAM desde Biopython.

Modularización realizada para aislar la lógica de alineamiento del resto del pipeline, permitiendo pruebas independientes y reutilización flexible.
"""
from Bio.Align import PairwiseAligner
from Bio.Align import substitution_matrices
import time
import itertools

def cargar_matriz_sustitucion(nombre: str) -> dict:
    """
    Carga una matriz de sustitución de Biopython.

    Parámetros:
        nombre (str): Nombre de la matriz ('blosum62', 'pam250', etc.)

    Retorna:
        tuple: matriz de sustitución.
    """
    try:
        return substitution_matrices.load(nombre)
    except Exception:
        raise ValueError(f"No se pudo cargar la matriz '{nombre}'.")


def smith_waterman(seq1: str, seq2: str, matrix_name: str, gap_penalty: int):
    """
    Implementación del algoritmo Smith-Waterman para alineamiento local.

    Retorna:
        Tuple[str, str, int]: secuencia1_alineada, secuencia2_alineada, puntuación
    """

    #Creamos la matriz inicial
    matriz = []
    for i in range(len(seq1) + 1):
        matriz.append([])
        for j in range(len(seq2) + 1):
            matriz[i].append(0)

    #Obtenemos la matriz de sustitucion, si la hay
    s_matrix = cargar_matriz_sustitucion(matrix_name)

    #Definimos la funcion score
    def score(a, b):
        try:
            return s_matrix[a, b]
        except KeyError:
            return s_matrix[b, a]
        except KeyError:
            return -1

    #Rellenamos la matriz
    for i in range(1,len(seq1) + 1):
        for j in range(1,len(seq2) + 1):
            matriz[i][j] = max(
                0,
                matriz[i - 1][j - 1] + score(seq1[i - 1], seq2[j - 1]),
                matriz[i - 1][j] - gap_penalty,
                matriz[i][j - 1] - gap_penalty
            )
    #Obtenemos la celda de mayor valor
    maximo = (0,0,0)
    for i in range(len(seq1) + 1):
        for j in range(len(seq2) + 1):
            if matriz[i][j] > maximo[0]:
                maximo = (matriz[i][j],i,j)
    # Empezamos desde la celda con el valor máximo
    score_final, i, j = maximo

    # Strings alineados (los construiremos al revés y luego invertimos)
    aln1 = ""
    aln2 = ""

    # Traceback mientras no llegamos a un 0
    while matriz[i][j] != 0:
        # Caracteres actuales de las secuencias
        letra1 = seq1[i - 1]
        letra2 = seq2[j - 1]

        # Puntuación por alinearlos
        puntaje_diag = matriz[i - 1][j - 1] + score(letra1, letra2)
        puntaje_arriba = matriz[i - 1][j] - gap_penalty
        puntaje_izquierda = matriz[i][j - 1] - gap_penalty

        # Caso 1: viene de la diagonal → match o mismatch
        if matriz[i][j] == puntaje_diag:
            aln1 = letra1 + aln1
            aln2 = letra2 + aln2
            i -= 1
            j -= 1

        # Caso 2: viene de arriba → gap en seq2
        elif matriz[i][j] == puntaje_arriba:
            aln1 = letra1 + aln1
            aln2 = "-" + aln2
            i -= 1

        # Caso 3: viene de izquierda → gap en seq1
        elif matriz[i][j] == puntaje_izquierda:
            aln1 = "-" + aln1
            aln2 = letra2 + aln2
            j -= 1

        else:
            # Seguridad: si no coincide con nada, rompemos (caso raro)
            break

    return aln1, aln2, score_final


def medir_rendimiento_custom(seqs: list, matrix_name: str, gap_penalty: int):
    """
    Mide el tiempo de ejecución del algoritmo propio Smith-Waterman.

    Parámetros:
        seqs (list): Lista de SeqRecord
    """
    start = time.perf_counter()
    for s1, s2 in itertools.combinations(seqs, 2):
        smith_waterman(str(s1.seq), str(s2.seq), matrix_name, gap_penalty)
    end = time.perf_counter()
    print(f"Tiempo ejecución Smith-Waterman propio: {end - start:.4f} segundos")


def medir_rendimiento_biopython(seqs: list, matrix_name: str):
    """
    Mide el tiempo de ejecución del alineador de Biopython.
    """
    # Configurar el alineador
    aligner = PairwiseAligner()
    aligner.mode = 'local'  # alineamiento local
    aligner.substitution_matrix = substitution_matrices.load(matrix_name)
    aligner.open_gap_score = -10
    aligner.extend_gap_score = -1

    # Medir el tiempo
    start = time.perf_counter()
    for s1, s2 in itertools.combinations(seqs, 2):
        alignments = aligner.align(str(s1.seq), str(s2.seq))
    end = time.perf_counter()
    print(f"Tiempo ejecución Smith-Waterman Biopython: {end - start:.4f} segundos")


def guardar_resultados_alineamientos(resultados: list, archivo_salida: str):
    """
    Guarda en un archivo los alineamientos obtenidos con Smith-Waterman propio.

    Cada resultado debe ser una tupla: (id1, id2, alineado1, alineado2, score)
    """
    with open(archivo_salida, "w") as f:
        for id1, id2, aln1, aln2, score in resultados:
            f.write(f">{id1} vs {id2} | Score: {score}\n")
            f.write(f"{aln1}\n{aln2}\n\n")