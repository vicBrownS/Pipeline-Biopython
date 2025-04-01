import os
from collections.abc import Iterator

from Bio.Seq import Seq
from Bio.SeqUtils import gc_fraction

import re
from itertools import zip_longest
from Bio import SeqIO
from Bio.SeqIO import SeqRecord
import logging

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

def guardar_records(lista_paths: List[str], lista_format: List[str], output_name: str):
    """
    Carga secuencias desde múltiples archivos y guarda todas las secuencias combinadas en un único archivo FASTA.

    Parámetros:\n
    - lista_paths (List[str]): Lista de rutas de los archivos que contienen las secuencias a cargar.
    - lista_format (List[str]): Lista de formatos correspondientes a cada archivo (por ejemplo, 'fasta', 'fastq', 'genbank').
      La longitud de esta lista debe coincidir con la de `lista_paths`.
    - output_name (str): Nombre del archivo de salida donde se guardarán las secuencias combinadas (sin extensión).

    Funcionamiento:\n
    1. Verifica que las listas proporcionadas no sean nulas, tengan la misma longitud y contengan solo cadenas.
    2. Si `output_name` no se proporciona, se asigna por defecto el nombre `"Secuencias"`.
    3. Itera sobre las listas de archivos y formatos:
       - Carga las secuencias usando `load_sequences()`.
       - Si la secuencia cargada es una instancia de `SeqRecord`, la añade directamente a la lista `secuencias`.
       - Si es un iterable (como una lista de `SeqRecord`), extiende la lista con todas las secuencias cargadas.
    4. Guarda todas las secuencias combinadas en un archivo FASTA utilizando `save_sequences()`.

    Excepciones:\n
    - ValueError: Si las listas son nulas, tienen longitudes diferentes o contienen elementos no válidos.
    - TypeError: Si algún elemento de las listas no es una cadena.

    Notas: \n
    - Esta función permite fusionar secuencias de diferentes archivos y formatos en un solo archivo.
    - Si se cargan secuencias duplicadas, se guardarán tal como se cargaron a menos que se eliminen previamente.
    - Se recomienda revisar las secuencias antes de guardarlas si se desea evitar duplicados.
    """

    #Verificaciones de validez de las listas de entrada
    if lista_paths is None or lista_format is None:
        raise ValueError("Las listas no pueden ser nulas")

    if not len(lista_format) == len(lista_paths):
        raise ValueError("Las listas deben ser de la misma longitud")

    if not all(isinstance(path, str) for path in lista_paths) or not all(
            isinstance(fmt, str) for fmt in lista_format):
        raise ValueError("Las listas solo deben contener cadenas")

    #Asignación de nombre de salida por defecto si no se proporciona
    if output_name is None:
        output_name = "Secuencias"

    # Lista para almacenar todas las secuencias combinadas
    secuencias = list()

    #Cargar secuencias desde cada archivo especificado
    for path,format in zip(lista_paths,lista_format):
        secuencia = load_sequences(path, format, verbose=False)

        #Si la carga devuelve una sola secuencia (SeqRecord), añádela directamente
        if isinstance(secuencia, SeqRecord):
            secuencias.append(secuencia)

        #Si devuelve múltiples secuencias (iterable), extiende la lista
        else:
            secuencias.extend(secuencia)

    #Guardar todas las secuencias cargadas en un único archivo FASTA
    save_sequences(secuencias, f"data/{output_name}")

    logging.info(f"Archivo {output_name}.fasta creado correctamente")

def convertir_formato(lista_secuencias,formato_salida: str, output_name: str):
    """
    Convierte las secuencias cargadas a un nuevo formato y las guarda en un archivo.

    Parámetros:
    - lista_secuencias (list): lista de secuencias a convertir
    - formato_salida (str): Formato de salida ('fasta', 'genbank', 'fastq', etc.).
    - output_name (str): Nombre del archivo de salida sin la extensión.

    Retorna:
    - None (guarda el archivo en el formato especificado).

    Excepciones:
    - ValueError: Si el formato de salida no es válido.
    - ValueError: Si `self.lista_secuencias` está vacía.
    """
    # Validación del formato de salida
    if not formato_salida or formato_salida.strip() == "":
        raise ValueError("El formato de salida debe tener algún valor.")

    # Mapeo de formatos compatibles con Biopython
    formatos_validos = {"fasta", "genbank", "fastq"}

    # Corregir formato en caso de alias (por ejemplo, "gb" → "genbank")
    formato_salida = "genbank" if formato_salida == "gb" else formato_salida.lower()

    if formato_salida not in formatos_validos:
        raise ValueError(f"El formato de salida '{formato_salida}' no es válido. Opciones: {formatos_validos}")

    # Verificación de lista vacía
    if not lista_secuencias:
        raise ValueError("No hay secuencias cargadas para convertir.")

    # Asegurar que el nombre del archivo tenga la extensión correcta
    extension_dict = {"fasta": ".fasta", "genbank": ".gb", "fastq": ".fastq"}
    extension = extension_dict[formato_salida]
    output_name = output_name.strip()
    if not output_name.endswith(extension):
        output_name += extension

    # Solución para GenBank: Añadir molecule_type si falta
    if formato_salida == "genbank":
        for seq in lista_secuencias:
            if "molecule_type" not in seq.annotations:
                if str(seq.seq).upper().find("U") == -1:
                    seq.annotations["molecule_type"] = "DNA"
                else:
                    seq.annotations["molecule_type"] = "RNA"

    # Solución para FASTQ: Añadir valores de calidad ficticios si faltan
    if formato_salida == "fastq":
        for seq in lista_secuencias:
            if "phred_quality" not in seq.letter_annotations:
                seq.letter_annotations["phred_quality"] = [40] * len(seq.seq)  # Calidad máxima ficticia

    # Guardar las secuencias en el formato especificado
    SeqIO.write(lista_secuencias, f"output/{output_name}", formato_salida)

    logging.info(f"Archivo guardado como '{output_name}' en formato '{formato_salida}'.")