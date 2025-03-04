import os
from collections.abc import Iterator

from Bio.Seq import Seq
from Bio.SeqUtils import gc_fraction
from Bio.Data import CodonTable
import re
from itertools import zip_longest
from Bio import SeqIO
from Bio.SeqIO import SeqRecord
import logging


def count_motifs(dna_seq: Seq, motifs: str or list, overlapping: bool) -> int or dict:
    """
    Cuenta la frecuencia de uno o varios patrones (motifs) en una secuencia de ADN.

    Parámetros:
    - dna_seq (Seq): Secuencia de ADN donde se busca el patrón.
    - motifs (str or list): Patrón de ADN (str) o lista de patrones a contar.
    - overlapping (bool): Indica si se cuentan ocurrencias con solapamiento.

    Retorna:
    - int o dict: Si es un solo motivo, retorna el número de apariciones.
                  Si es una lista, retorna un diccionario con las frecuencias de cada motivo.

    Excepciones:
    - ValueError: Si la secuencia o los motivos son nulos.
    """

    # Verifica que la secuencia de ADN no sea nula
    if dna_seq is None:
        raise ValueError("La secuencia de ADN no puede ser nula")

    # Verifica que se haya proporcionado al menos un motivo
    if motifs is None:
        raise ValueError("Inserte algún motivo")

    # Si los motivos son una lista, contar cada uno y devolver un diccionario
    if isinstance(motifs, list):
        resultados = {}  # Diccionario para almacenar los conteos de cada motivo
        for motif in motifs:
            if overlapping:
                # Cuenta las ocurrencias con solapamiento y lo almacena en el diccionario
                resultados[motif] = dna_seq.count_overlap(motif.upper())
            else:
                # Cuenta las ocurrencias sin solapamiento y lo almacena en el diccionario
                resultados[motif] = dna_seq.count(motif.upper())
        return resultados  # Devuelve el diccionario con los conteos

    else:
        # Si es un solo motivo, se cuenta directamente
        if overlapping:
            return dna_seq.count_overlap(motifs.upper())
        else:
            return dna_seq.count(motifs.upper())


def gc_content(dna_seq: Seq, pos_start: int = 0, pos_end: int = None) -> float:
    """
    Calcula el porcentaje de nucleótidos G y C en una secuencia de ADN,
    ya sea en toda la secuencia o en un rango específico.

    Parámetros:
    - dna_seq (Seq): Secuencia de ADN.
    - pos_start (int, opcional): Posición de inicio (por defecto es 0).
    - pos_end (int, opcional): Posición de fin (por defecto es el final de la secuencia).

    Retorna:
    - float: Porcentaje de nucleótidos G y C en el rango especificado.

    Excepciones:
    - ValueError: Si la secuencia es nula.
    - ValueError: Si las posiciones de inicio o fin son inválidas.
    """

    # Verifica que la secuencia de ADN no sea nula
    if dna_seq is None:
        raise ValueError("La secuencia de ADN no puede ser nula")

    # Si no se proporciona una posición de fin, se asigna la longitud de la secuencia
    if pos_end is None:
        pos_end = len(dna_seq)

    # Si se especifica la posicion de comienzo como nula, se le asigna la primera posición.
    if pos_start is None:
        pos_start = 0

    # Verificación de errores en las posiciones de inicio y fin
    if (
        pos_start > pos_end or  # El inicio no puede ser mayor que el final
        pos_start < 0 or  # No se permiten índices negativos
        pos_end < 0 or  # No se permite que el fin sea negativo
        pos_end > len(dna_seq)  # No se permite que el fin sea mayor que la secuencia
    ):
        raise ValueError("Error en la posición de comienzo o posición final")

    # Calcula la fracción de nucleótidos G y C en el rango especificado
    fraccion_gc = gc_fraction(dna_seq[pos_start:pos_end])

    # Retorna el porcentaje redondeado a dos decimales
    return round(fraccion_gc * 100, 2)

def get_coding_strand(mRNA_seq: Seq) -> Seq:
    """
    Obtiene la hebra codificante de ADN a partir de una secuencia de ARN mensajero (mRNA).

    Parámetros:
    - mRNA_seq (Seq): Secuencia de ARN mensajero.

    Retorna:
    - Seq: Hebra codificante de ADN (idéntica al mRNA, con 'U' reemplazado por 'T').

    Excepciones:
    - ValueError: Si la secuencia es nula o contiene nucleótidos no válidos.
    """

    # Verificación de que la secuencia no sea nula
    if mRNA_seq is None:
        raise ValueError("La secuencia de ARN no puede ser nula.")

    # Verificación de nucleótidos válidos
    if not set(str(mRNA_seq.upper())).issubset({"A", "U", "G", "C"}):
        raise ValueError("La secuencia de ARN contiene nucleótidos no válidos.")

    # Reemplaza uracilos (U) por timinas (T) para obtener la hebra codificante de ADN
    dna_codificante = mRNA_seq.replace("U", "T")

    return dna_codificante

def identify_codons(mRNA_Seq: Seq, table: str or int = "Standard"):
    """
    Identifica los codones de inicio y parada en una secuencia de ARN mensajero (mRNA),
    según la tabla genética especificada, y devuelve las posibles secuencias codificantes.

    Parámetros:
    - mRNA_Seq (Seq): Secuencia de ARN mensajero.
    - table (str or int): Nombre o ID de la tabla genética (por defecto: "Standard").

    Retorna:
    - list: Lista de tuplas con ((inicio, fin), secuencia codificante).

    Excepciones:
    - ValueError: Si la secuencia o tabla son nulas, o si la secuencia contiene nucleótidos no válidos.
    - KeyError: Si la tabla genética proporcionada no es válida.
    """

    # Verificación de que la secuencia y la tabla no sean nulas
    if table is None or mRNA_Seq is None:
        raise ValueError("La secuencia de ARN no puede ser nula, ni tampoco la tabla")

    # Validación de que la secuencia solo contiene nucleótidos válidos
    if not set(str(mRNA_Seq.upper())).issubset({"A", "U", "G", "C"}):
        raise ValueError("La secuencia de ARN contiene nucleótidos no válidos")

    # Obtener la tabla genética proporcionada (por nombre o por ID)
    try:
        if isinstance(table, str):
            codones = CodonTable.unambiguous_rna_by_name[table]  # Corrección: Se debe usar 'unambiguous_rna_by_name'
        elif isinstance(table, int):
            codones = CodonTable.unambiguous_rna_by_id[table]    # Corrección: Se debe usar 'unambiguous_rna_by_id'
    except KeyError:
        raise KeyError("Nombre o ID de tabla incorrecto")

    # Obtener codones de inicio y parada de la tabla genética
    codones_inicio = codones.start_codons
    codones_parada = codones.stop_codons

    # Buscar todas las posiciones de inicio y parada en la secuencia
    posiciones_inicio = [match.start() for cod_in in codones_inicio for match in re.finditer(cod_in, str(mRNA_Seq))]
    posiciones_parada = [match.end() for cod_pa in codones_parada for match in re.finditer(cod_pa, str(mRNA_Seq))]

    # Extraer secuencias codificantes entre codones de inicio y parada
    secuencias = []
    for posicion_inicio in posiciones_inicio:
        for posicion_parada in posiciones_parada:
            if posicion_inicio + 2 < posicion_parada - 2: # Asegura que la relación entre las posciones de incio y parada sea la correcta
                secuencia = mRNA_Seq[posicion_inicio:posicion_parada]
                if len(secuencia) % 3 == 0: #Comprueba que la secuencia este correctamente formada
                    secuencias.append(((posicion_inicio, posicion_parada), mRNA_Seq[posicion_inicio:posicion_parada]))

    return secuencias


def check_translation(dna: Seq, peptide: Seq, table: int or str = "Standard", cds: list = None,
                      to_stop: bool = True) -> bool or list:
    """
    Verifica si la traducción de una secuencia de ADN coincide con la secuencia de aminoácidos proporcionada.

    Parámetros:
    - dna (Seq): Secuencia de ADN a traducir.
    - peptide (Seq): Secuencia de aminoácidos esperada (proporcionada por GenBank).
    - table (str or int): Nombre o ID de la tabla genética utilizada para la traducción (por defecto: "Standard").
    - cds (list, opcional): Lista que indica el rango de la región codificante (inicio, fin). Si no se proporciona, se traduce toda la secuencia.
    - to_stop (bool): Si es True, la traducción se detendrá en el primer codón de parada (por defecto: True).

    Retorna:
    - bool: True si ambas secuencias coinciden completamente.
    - list: Lista de tuplas con las posiciones y los aminoácidos que no coinciden si hay diferencias.
      Cada tupla tiene la forma (posición, aminoácido traducido, aminoácido esperado).

    Excepciones:
    - ValueError: Si alguna de las secuencias, la tabla o el rango de CDS son nulos o inválidos.
    - KeyError: Si el nombre o ID de la tabla genética no es válido.
    """

    # Verificación de entradas nulas
    if dna is None or peptide is None or table is None:
        raise ValueError("Secuencia de ADN, péptido o tabla genética no pueden ser nulas.")

    # Asignación por defecto de CDS si no se proporciona
    if cds is None:
        cds = [0, len(dna)]

    # Verificacion de que el rango cds tiene longitud y valores correctos
    if not isinstance(cds, list) or len(cds) != 2 or not all(isinstance(i, int) for i in cds):
        raise ValueError("El parámetro cds debe ser una lista de dos enteros.")

    # Verificación de que el rango CDS sea válido
    if cds[0] < 0 or cds[1] > len(dna) or cds[0] > cds[1]:
        raise ValueError("Valores de rango de CDS incorrectos.")

    # Obtener la tabla genética proporcionada (por nombre o por ID)
    try:
        if isinstance(table, str):
            codones = CodonTable.unambiguous_dna_by_name[table]  # Accede a la tabla genética por nombre
        elif isinstance(table, int):
            codones = CodonTable.unambiguous_dna_by_id[table]  # Accede a la tabla genética por ID
    except KeyError:
        raise KeyError("Nombre o ID de tabla genética incorrecto.")

    # Traducción del segmento de ADN especificado por CDS hasta el primer codón de parada (si to_stop=True)
    supposed_peptide = dna[cds[0]:cds[1]].translate(table=codones, to_stop=to_stop)

    # Si las secuencias coinciden, devuelve True directamente
    if peptide == supposed_peptide:
        return True
    else:
        # Si no coinciden, se generan las posiciones de las diferencias
        lista_no_coincidencias = [
            (i, aa_traducido or "-", aa_esperado or "-")  # Maneja valores faltantes con "-"
            for i, (aa_traducido, aa_esperado) in enumerate(zip_longest(supposed_peptide, peptide))
            if aa_traducido != aa_esperado
        ]
        return lista_no_coincidencias

def load_sequences(file_path: str, file_format: str = None, verbose: bool = False) -> list:
    """
    Carga y explora secuencias desde un archivo en formatos FASTA, FASTQ o GenBank.

    Parámetros:
    - file_path (str): Ruta al archivo que contiene las secuencias.
    - file_format (str, opcional): Formato del archivo ('fasta', 'fastq', 'gb'). Si no se especifica, se infiere de la extensión del archivo.
    - verbose (bool, opcional): Si es True, imprime un resumen detallado de las secuencias cargadas.

    Retorna:
    - SeqRecord: Un iterador de objetos SeqRecord con las secuencias cargadas.

    Excepciones:
    - ValueError: Si el formato es inválido o si la ruta del archivo es None.
    """

    # Validación de que la ruta del archivo no sea None
    if file_path is None:
        raise ValueError("El camino del archivo no puede ser None")

    # Si no se especifica el formato, se infiere de la extensión del archivo
    if file_format is None:
        file_format = file_path.split(".")[-1]  # Obtiene la extensión del archivo
        if file_format == "fas":  # Maneja la extensión .fas como sinónimo de .fasta
            file_format = "fasta"

    #Verificación de que el formato sea válido
    formatos_validos = {"fasta", "fastq", "gb"}
    if file_format not in formatos_validos:
        raise ValueError("Formato inválido. Los formatos aceptados son: fasta, fastq, gb")

    # Carga las secuencias con SeqIO
    try:
        seq_record = SeqIO.parse(file_path, file_format)
        if not seq_record:
            raise ValueError("El archivo no contiene secuencias o está vacío.")
    except FileNotFoundError:
        raise FileNotFoundError(f"Archivo no encontrado: {file_path}")
    # Convertir a lista para iterar más de una vez si es necesario
    secuencias = list(seq_record)

    # Si verbose es True, imprime información detallada sobre las secuencias
    # Dentro de la función donde se encuentra el bloque original:
    if verbose:
        logging.info(f"\nNúmero total de secuencias: {len(secuencias)}\n")

        for i, secuencia in enumerate(secuencias, 1):
            logging.info(f"Secuencia {i}:")
            logging.info(f"  ID: {secuencia.id}")
            logging.info(f"  Nombre: {secuencia.name}")
            logging.info(f"  Descripción: {secuencia.description}")
            logging.info(f"  Longitud: {len(secuencia.seq)} nucleótidos")
            logging.info(
                f"  Secuencia (primeros 50 caracteres): {secuencia.seq[:50]}{'...' if len(secuencia.seq) > 50 else ''}")

            if file_format in ["fastq", "gb"]:
                logging.info(f"  Anotaciones: {secuencia.annotations or 'Vacío'}")
                logging.info(f"  Anotaciones de letras: {secuencia.letter_annotations or 'Vacío'}")

            logging.info("-" * 50)

    return secuencias

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

def extract_subsequences(sequences: Iterator[SeqRecord] or list, start_motif: str, end_motif: str) -> dict:
    """
    Extrae subsecuencias comprendidas entre los motivos de inicio y fin, incluyendo dichos motivos.

    Parámetros:
    - sequences (Iterator[SeqRecord] or list): Conjunto de secuencias a procesar.
    - start_motif (str): Motivo de inicio para la extracción.
    - end_motif (str): Motivo de finalización para la extracción.

    Retorna:
    - dict: Diccionario con identificadores de secuencias como claves y las subsecuencias extraídas como valores.

    Excepciones:
    - ValueError: Si las secuencias o motivos son nulos o inválidos.
    - ValueError: Si las secuencias no son instancias de SeqRecord.
    """

    secuencias = dict()

    #Verificación de que los motivos y las secuencias no sean nulos
    if start_motif is None or end_motif is None or sequences is None:
        raise ValueError("Valores de motivos o secuencias nulos")

    #Validacion de motivos
    valid_nucleotides = {"A", "C", "G", "T", "U"}
    if not set(start_motif.upper()).issubset(valid_nucleotides) or not set(end_motif.upper()).issubset(valid_nucleotides):
        raise ValueError("Motivos de comienzo o fin no válidos. Deben contener solo A, C, G, T, o U.")

    #Verificación de que todas las secuencias sean instancias de SeqRecord
    if not all(isinstance(seq, SeqRecord) for seq in sequences):
        raise ValueError("Todas las secuencias tienen que ser instancias de SeqRecord.")

    #Búsqueda y extracción de subsecuencias
    for seq in sequences:
        # Encuentra la primera aparición del motivo de inicio
        start_pos = str(seq.seq).find(start_motif)
        # Encuentra la primera aparición del motivo de finalización después del motivo de inicio
        end_pos = str(seq.seq).find(end_motif, start_pos)

        if start_pos != -1 and end_pos != -1 and end_pos >= start_pos:
            subsecuencia = seq.seq[start_pos:end_pos + len(end_motif)]
            secuencias[seq.id] = subsecuencia

    return secuencias

def index_sequences(file_path: str, file_format: str = None) -> dict:
    """
    Crea un diccionario indexado de secuencias basado en sus identificadores utilizando SeqIO.index().

    Parámetros:
    - file_path (str): Ruta al archivo de secuencias.
    - file_format (str, opcional): Formato del archivo ('fasta', 'fastq', 'genbank', 'gb').
      Si no se especifica, se infiere de la extensión.

    Retorna:
    - Mapping[str, SeqIO.SeqRecord]: Diccionario indexado con IDs como claves y SeqRecords como valores.

    Excepciones:
    - ValueError: Si la ruta es nula o el formato es inválido.
    - FileNotFoundError: Si el archivo no existe.
    - Exception: Para otros errores durante el indexado.
    """

    #Verificación de ruta
    if file_path is None:
        raise ValueError("La ruta del archivo no puede ser nula.")
    if not os.path.isfile(file_path):
        raise FileNotFoundError(f"Archivo no encontrado: {file_path}")

    #Inferencia y validación del formato
    if file_format is None:
        file_format = file_path.split(".")[-1]
        if file_format == "fas":
            file_format = "fasta"

    formatos_validos = {"fasta", "fastq", "genbank", "gb"}
    if file_format not in formatos_validos:
        raise ValueError(f"Formato '{file_format}' no válido. Formatos aceptados: {formatos_validos}")

    #Indexación
    try:
        return SeqIO.index(file_path, file_format)
    except Exception as e:
        raise Exception(f"Error al indexar el archivo: {e}")

def get_sequences_by_id(seq_id: str, seqs_index) -> SeqRecord:
    """
    Recupera una secuencia a partir de un identificador utilizando un índice de secuencias.

    Parámetros:
    - seq_id (str): Identificador de la secuencia a recuperar.
    - seqs_index (Mapping[str, SeqRecord]): Diccionario indexado de secuencias.

    Retorna:
    - SeqRecord: Secuencia correspondiente al identificador proporcionado.

    Excepciones:
    - ValueError: Si el identificador es nulo, vacío o el índice es nulo/vacío.
    - KeyError: Si el identificador no se encuentra en el índice.
    """

    #Verificación de entradas
    if not isinstance(seq_id, str) or not seq_id.strip():
        raise ValueError("El identificador debe ser una cadena no vacía.")
    if len(seqs_index) == 0:
        raise ValueError("El índice de secuencias está vacío.")

    #Acceso a la secuencia
    try:
        return seqs_index[seq_id]
    except KeyError:
        raise KeyError(f"ID '{seq_id}' no encontrado. IDs disponibles: {list(seqs_index.keys())}")

def save_sequences(sequences: list or Iterator[SeqRecord], output_file_name: str):
    """
    Guarda un conjunto de secuencias en un archivo FASTA y muestra estadísticas.

    Parámetros:
    - sequences (list or Iterator[SeqRecord]): Lista o iterador de secuencias.
    - output_file_name (str): Nombre del archivo de salida (sin extensión).

    Funcionalidades:
    - Guarda las secuencias en formato FASTA.
    - Muestra el número de secuencias almacenadas y la longitud media.

    Excepciones:
    - ValueError: Si las secuencias son nulas, inválidas o vacías.
    """

    #Verificaciones
    if sequences is None:
        raise ValueError("Las secuencias no pueden ser nulas.")
    if not isinstance(sequences, (list, Iterator)):
        raise ValueError("Entrada de secuencias inválida. Debe ser una lista o un iterador.")
    if len(sequences) == 0:
        raise ValueError("La lista de secuencias está vacía.")
    if not all(isinstance(seq, SeqRecord) for seq in sequences):
        raise ValueError("Todas las secuencias deben ser instancias de SeqRecord.")

    #Asignación de nombre de archivo
    if not output_file_name or not output_file_name.strip():
        output_file_name = "Secuencias"
    if not output_file_name.endswith(".fasta"):
        output_file_name += ".fasta"

    #Escritura de las secuencias
    SeqIO.write(sequences, output_file_name, "fasta")

    #Estadísticas
    num_secuencias = len(sequences)
    longitud_media = sum(len(seq.seq) for seq in sequences) / num_secuencias

    #Mensajes informativos
    print(f"Secuencias guardadas en '{output_file_name}'")
    print(f"Número de secuencias almacenadas: {num_secuencias}")
    print(f"Longitud media de las secuencias: {round(longitud_media,3)}")





