from typing import List, Dict, Iterator
from Bio.SeqUtils import nt_search
from Bio.SeqIO import SeqRecord
from Bio.Seq import Seq



def buscar_motivos(self, motivos: List[str]) -> Dict[str, Dict[str, List[int]]]:
    """
    Busca la posición de uno o más motivos en cada secuencia.

    Parámetros:
    - motivos (List[str]): Lista de secuencias de ADN a buscar.

    Retorna:
    - Dict[str, Dict[str, List[int]]]:
      Diccionario con los IDs de las secuencias y otro diccionario con:
        - Claves: Motivos buscados.
        - Valores: Lista con las posiciones donde aparece el motivo.

    Excepciones:
    - ValueError: Si la lista de motivos está vacía o contiene elementos no válidos.
    """

    # Validaciones de entrada
    if not motivos:
        raise ValueError("La lista de motivos no puede estar vacía.")

    if not all(isinstance(motivo, str) for motivo in motivos):
        raise ValueError("Todos los motivos deben ser cadenas de texto (str).")

    # Diccionario para almacenar los resultados
    motivos_dict = {}

    # Buscar motivos en cada secuencia
    for record in self.lista_secuencias:
        motivos_dict[record.id] = {}
        for motivo in motivos:
            posiciones = nt_search(str(record.seq), motivo.upper())
            # Si el motivo se encuentra, almacenar posiciones
            if len(posiciones) > 1:
                motivos_dict[record.id][motivo] = posiciones[1:]
            else:
                motivos_dict[record.id][motivo] = []  # Si no se encuentra, lista vacía

    return motivos_dict

def obtener_estadisticas(self) -> dict:
    """
    Calcula estadísticas generales sobre las secuencias en `lista_secuencias`.

    Retorna:
    - dict con:
      - 'total_secuencias': Número total de secuencias.
      - 'longitud_max': Longitud de la secuencia más larga.
      - 'longitud_min': Longitud de la secuencia más corta.
      - 'longitud_media': Promedio de longitud de las secuencias.
      - 'contenido_gc_medio': Promedio del porcentaje de GC en todas las secuencias.
    """

    total_secuencias = len(self.lista_secuencias)

    # Manejo de lista vacía
    if total_secuencias == 0:
        return {"total_secuencias": 0, "longitud_max": None, "longitud_min": None,
                "longitud_media": None, "contenido_gc_medio": None}

    # Recorrer la lista solo una vez para calcular todas las estadísticas
    total_longitudes = 0
    total_gc = 0
    longitud_max = -1.0
    longitud_min = 100000.0

    for secuencia in self.lista_secuencias:
        len_seq = len(secuencia.seq)
        total_longitudes += len_seq
        total_gc += gc_content(secuencia.seq)

        if len_seq > longitud_max:
            longitud_max = len_seq
        if len_seq < longitud_min:
            longitud_min = len_seq

    #Calcular valores finales
    longitud_media = total_longitudes / total_secuencias
    contenido_gc_medio = total_gc / total_secuencias

    return {"total_secuencias": total_secuencias, "longitud_max": longitud_max, "longitud_min": longitud_min,
            "longitud_media": longitud_media, "contenido_gc_medio": contenido_gc_medio}

def filter_sequences_by_length(min_length: int, max_length: int, sequences: Iterator[SeqRecord] or list) -> Iterator[SeqRecord]:
    """
    Filtra un conjunto de secuencias y retorna solo aquellas que estén dentro del rango de longitud especificado.

    Parámetros:
    - min_length (int): Longitud mínima permitida.
    - max_length (int): Longitud máxima permitida.
    - sequences (Iterator[SeqRecord] or list): Conjunto de secuencias a filtrar.

    Retorna:
    - Iterator[SeqRecord]: Iterador con las secuencias filtradas.

    Excepciones:
    - ValueError: Si las secuencias son nulas o las longitudes son inválidas.
    - TypeError: Si las secuencias no son instancias de SeqRecord.
    """

    #Validaciones básicas
    if sequences is None:
        raise ValueError("Las secuencias no deben ser nulas.")
    if not all(isinstance(seq, SeqRecord) for seq in sequences):
        raise TypeError("Todas las secuencias deben ser instancias de SeqRecord.")

    if min_length is None:
        min_length = 0
    if max_length is None:
        max_length = max(len(seq) for seq in sequences)  # Longitud máxima de las secuencias

    if min_length >= max_length or min_length < 0:
        raise ValueError("Longitudes de secuencia incorrectas.")

    #Retorno directo de un generador
    return (seq for seq in sequences if min_length <= len(seq) <= max_length)

