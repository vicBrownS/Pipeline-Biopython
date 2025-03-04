from Utils import *
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from collections import Counter
from typing import List, Dict, Tuple
from Bio.SeqUtils import nt_search
from Bio.Data import CodonTable
import logging


class Pipeline:
    """
    Clase que gestiona el procesamiento de secuencias biológicas (ADN, ARN y proteínas).

    **Funcionalidades principales:** \n
    1. **Carga de Secuencias:** \n
       - `load_sequences(file_path, file_format)`: Carga secuencias desde archivos FASTA, GenBank o FASTQ.
       - `add_sequences(sequences, reindex=True)`: Añade secuencias a la lista manualmente.

    2️. **Manejo de Secuencias:** \n
       - `get_sequence_by_id(seq_id)`: Obtiene una secuencia específica por su ID.
       - `get_sequence_by_name(name)`: Obtiene una secuencia por su nombre.
       - `eliminar_secuencias_by_id(seq_id, all=False)`: Elimina secuencias específicas o todas a la vez.
       - `eliminar_duplicados(id=True, name=True)`: Elimina secuencias duplicadas por ID o nombre.
       - `index_sequences_by_id_name(force=False)`: Crea índices para acceso rápido a las secuencias.

    3. **Procesamiento de Secuencias:** \n
       - `filtrar_secuencias_por_longitud(min_len, max_len, reindex=True)`: Filtra secuencias según su longitud.
       - `buscar_motivos(motivos)`: Busca posiciones de motivos dentro de las secuencias.
       - `obtener_secuencia_sin_exones(seq_id)`: Extrae solo la región codificante (sin intrones).
       - `asignar_cds(seq_id, inicio, final)`: Permite asignar manualmente una región CDS.

    4. **Análisis de Secuencias:** \n
       - `obtener_estadisticas()`: Calcula estadísticas como longitud media y contenido GC.

    5. **Traducción de Secuencias:** \n
       - `traducir_secuencias(table="Standard", to_stop=False)`: Traduce secuencias de ADN a proteínas.

    6. **Conversión y Guardado:** \n
       - `convertir_formato(formato_salida, output_name)`: Convierte secuencias a FASTA, GenBank o FASTQ.
       - `guardar_records(lista_paths, lista_format, output_name)`: Carga múltiples archivos y los guarda en un único FASTA.

    **Notas:** \n
        - Utiliza `SeqRecord` de Biopython para almacenar y manipular secuencias.
        - Permite cargar secuencias desde múltiples archivos y combinarlas en un solo archivo FASTA.
        - Soporta traducción con diferentes tablas genéticas y extracción de regiones codificantes (CDS, exones).
        - Indexa secuencias para acceso rápido y permite filtrado por longitud o búsqueda de motivos.
    """


    def __init__(self):
        """
        Inicializa una nueva instancia de la clase Pipeline.

        Atributos:\n
        - lista_secuencias (List[SeqRecord]): Lista que almacena todas las secuencias cargadas o añadidas.
        - index_by_id (dict): Diccionario que permite acceder a las secuencias por su identificador (ID).
        - index_by_name (dict): Diccionario que permite acceder a las secuencias por su nombre.

        Notas:\n
        - Los índices se pueden actualizar llamando a `index_sequences_by_id_name()`.
        - `lista_secuencias` puede contener duplicados, pero los índices no.
        """
        self.lista_secuencias: List[SeqRecord] = []  # Lista de SeqRecords
        self.index_by_id = dict() # Índice por ID
        self.index_by_name = dict() # Índice por nombre
        self.cds_dict = dict() # {ID_secuencia: (inicio, fin)}
        self.exones_dict= dict() # {ID_secuencia: [(inicio, fin), (inicio2, fin2)]}

        # Configurar el logger
        logging.basicConfig(
            filename="pipeline.log",
            level=logging.INFO,
            format="%(asctime)s - %(levelname)s - %(message)s",
            datefmt="%Y-%m-%d %H:%M:%S",
            encoding = "utf-8"
        )

    def __str__(self):
        """
        Representación de las secuencias en la clase Pipeline.
        Retorna: \n
        - str que representa las secuencias que contiene la clase
        """
        string = ""
        string += "\n Secuencias:\n"
        for seq in self.lista_secuencias:
            string += f"ID: {seq.id} | Nombre: {seq.name} | Longitud: {len(seq.seq)} | Secuencia: {seq.seq[0:50]}...\n"
        return string

    def get_sequence_by_id(self, id: str):
        """
        Recupera una secuencia específica a partir de su identificador (ID).

        Parámetros:\n
        - id (str): El identificador de la secuencia que se desea recuperar.

        Retorna:\n
        - Seq: La secuencia correspondiente al identificador proporcionado.

        Excepciones:\n
        - ValueError: Se lanza si el identificador es nulo o inválido.
        - KeyError: Se lanza si el identificador no se encuentra en el índice de secuencias.

        Notas:\n
        - Se requiere que las secuencias estén indexadas previamente mediante `index_sequences_by_id_name()`.
        - Si el identificador no existe, se mostrará un mensaje de error.
        """
        try:
            seq = get_sequences_by_id(id, self.index_by_id)
            return seq
        except ValueError as v:
            logging.error(f"ValueError: {v}")
        except KeyError as k:
            logging.error(f"KeyError: {k}")

    def get_sequence_by_name(self, name: str):
        """
        Recupera una secuencia específica a partir de su nombre.

        Parámetros:\n
        - name (str): El nombre de la secuencia que se desea recuperar.

        Retorna:\n
        - Seq: La secuencia correspondiente al nombre proporcionado.

        Excepciones:\n
        - ValueError: Se lanza si el nombre es nulo o inválido.
        - KeyError: Se lanza si el nombre no se encuentra en el índice de secuencias.

        Notas:\n
        - Se requiere que las secuencias estén indexadas previamente mediante `index_sequences_by_id_name()`.
        - Si el nombre no existe, se mostrará un mensaje de error.
        """
        try:
            seq = get_sequences_by_id(name, self.index_by_name)
            return seq
        except ValueError as v:
            logging.error(f"ValueError: {v}")
        except KeyError as k:
            logging.error(f"KeyError: {k}")

    def add_sequences(self, sequences: list or Iterator[SeqRecord] or SeqRecord, reindex: bool = True):
        """
        Añade nuevas secuencias a la lista de secuencias del pipeline.

        Parámetros:\n
        - sequences (list o Iterator[SeqRecord]): Una lista o iterador de objetos `SeqRecord`
          que se desean agregar al conjunto de secuencias.
        - reindex (bool, opcional): Si se establece en `True`, se reindexarán automáticamente
          los diccionarios `index_by_id` y `index_by_name` después de añadir las secuencias.
          Por defecto es `False`.

        Retorna:\n
        - None

        Excepciones:\n
        - ValueError: Si las secuencias proporcionadas son nulas, vacías o no son del tipo correcto.
        - TypeError: Si algún elemento de la lista no es una instancia de `SeqRecord`.

        Notas:\n
        - Si `reindex` es `True`, se actualizarán los índices para que las nuevas secuencias
          estén disponibles para búsqueda inmediata.
        - Se pueden agregar secuencias duplicadas; utiliza `eliminar_duplicados()` si necesitas filtrarlas.
        """
        if sequences is None:
            raise ValueError("Las secuencias no pueden ser nulas")
        if not isinstance(sequences, (list, Iterator,SeqRecord)):
            raise ValueError("La entrada debe ser una lista o un iterador de secuencias")
        if len(sequences) == 0:
            raise ValueError("La lista de secuencias está vacía")
        if not all(isinstance(seq, SeqRecord) for seq in sequences) and not isinstance(sequences, SeqRecord):
            raise TypeError("Todas las secuencias deben ser instancias de SeqRecord")

        if isinstance(sequences, SeqRecord):
            #Añadimos secuencias individuales
            self.lista_secuencias.append(sequences)
        else:
            # Añadir las nuevas secuencias a la lista
            self.lista_secuencias.extend(sequences)


        # Reindexar si el parámetro reindex está activado
        if reindex:
            self.index_sequences_by_id_name()

    def load_sequences(self, file_path: str, file_format: str, verbose=True):
        """
        Carga secuencias de un archivo y extrae información sobre CDS y exones si el archivo es GenBank.

        Parámetros:
        - file_path (str): Ruta del archivo a cargar.
        - file_format (str): Formato del archivo ('fasta', 'genbank', 'fastq').
        - verbose (bool, opcional): Si es `True`, muestra información sobre las secuencias cargadas.

        Funcionamiento:
        1. Carga las secuencias en `lista_secuencias`.
        2. Si el formato es GenBank, extrae información de CDS y exones.
        3. Indexa las secuencias por ID y nombre.

        Excepciones:
        - ValueError: Si `file_path` o `file_format` son inválidos.
        - FileNotFoundError: Si el archivo no se encuentra.
        """
        try:
            # Carga secuencias desde el archivo y las añade sin reindexar individualmente
            for seq in load_sequences(file_path, file_format, verbose=verbose):
                self.add_sequences(seq, reindex=False)

                # Extraer información del CDS si el formato es GenBank
                if file_format.lower() == 'genbank' or file_format.lower() == 'gb':
                    for feature in seq.features:
                        if feature.type == 'CDS':
                            start, end = int(feature.location.start), int(feature.location.end)
                            self.cds_dict[seq.id] = (start, end)
                        elif feature.type == 'exon':
                            if seq.id not in self.exones_dict.keys():
                                self.exones_dict[seq.id] = []
                            self.exones_dict[seq.id].append((int(feature.location.start), int(feature.location.end)))

            # Reindexa después de cargar todas las secuencias para un acceso rápido
            self.index_sequences_by_id_name()

        except Exception as e:
            # Captura y muestra cualquier error ocurrido durante la carga o la indexación
            logging.error(f"Excepción: {e}")

    def get_sequences(self) -> list[SeqRecord]:
        """
        Devuelve la lista de secuencias que contiene el pipeline

        Retorna:\n
        - list: lista de secuencias (SeqRecord)
        """
        return self.lista_secuencias

    def get_sequences_id(self):
        """
        Devuelve las claves del diccionario

        Retorna:\n
        - list: claves de index_by_id
        """
        return self.index_by_id.keys()

    def get_sequences_name(self):
        return self.index_by_name.keys()

    def index_sequences_by_id_name(self, force=False):
        """
        Crea índices para acceder rápidamente a las secuencias por sus identificadores (`id`) y nombres (`name`).

        Parámetros:\n
        - force (bool, opcional): Si se establece en `True`, la indexación se realiza directamente,
          incluso si se detectan duplicados. Si es `False`, se pedirá confirmación al usuario antes de continuar.

        Funcionamiento:\n
        1. Busca duplicados de IDs y nombres en la lista de secuencias (`self.lista_secuencias`).
        2. Si se encuentran duplicados:
           - Se muestra una advertencia indicando los duplicados encontrados.
           - Se informa que en la indexación solo se conservará la última ocurrencia de cada duplicado.
           - Se aclara que los duplicados permanecerán en la lista original y se pueden eliminar más adelante.
           - Si `force` es `False`, se solicita confirmación al usuario para continuar.
        3. Si el usuario confirma o `force` es `True`, se crea el índice de secuencias por `id` y `name`.
           - En caso de duplicados, se conservará la última instancia encontrada.
        """

        # Busca duplicados de IDs en la lista de secuencias
        duplicates_id = [k for k, v in Counter(record.id for record in self.lista_secuencias).items() if v > 1]

        # Busca duplicados de nombres en la lista de secuencias
        duplicates_name = [k for k, v in Counter(record.name for record in self.lista_secuencias).items() if v > 1]

        # Flags para controlar si se debe proceder con la indexación
        flag_id = False if len(duplicates_id) > 0 else True  # False si hay duplicados de ID
        flag_name = False if len(duplicates_name) > 0 else True  # False si hay duplicados de nombre

        # Si se encuentran duplicados, se muestra la advertencia
        if not flag_id or not flag_name:
            logging.info(
                f"Advertencia. Duplicados encontrados de nombre o id.\n"
                f"Duplicados de id: {duplicates_id}\n"
                f"Duplicados de nombre: {duplicates_name}\n"
                f"Si hay duplicados solo se guardará la última instancia en los diccionarios indexados.\n"
                f"En la lista de secuencias se guardarán los registros duplicados.\n"
                f"Estos se podrán eliminar posteriormente."
            )

            # Si `force` es False, se solicita confirmación al usuario
            if not force:
                ad = input("Seguro que quiere continuar con la indexación (Y/N): ")
                if ad.lower() == "y":
                    flag_id = True
                    flag_name = True  # Se permite la indexación
                elif ad.lower() == "n":
                    logging.info("Indexación cancelada")
                    return  # Sale del método sin realizar la indexación

        # Si se confirma la indexación o no hay duplicados
        if flag_id and flag_name:
            # Crea el diccionario indexado por ID (última instancia con ID repetido prevalece)
            self.index_by_id = {seq.id: seq for seq in self.lista_secuencias}

            # Crea el diccionario indexado por nombre (última instancia con nombre repetido prevalece)
            self.index_by_name = {seq.name: seq for seq in self.lista_secuencias}

    def eliminar_duplicados(self, id: bool = True, name: bool = True):
        """
        Elimina secuencias duplicadas de la lista `lista_secuencias` según el criterio especificado.

        Parámetros:\n
        - id (bool, opcional): Si se establece en `True` (por defecto), se eliminarán duplicados basados en el identificador (`id`).
        - name (bool, opcional): Si se establece en `True` (por defecto), se eliminarán duplicados basados en el nombre (`name`).

        Retorna:\n
        - None

        Notas:\n
        - Si ambos parámetros (`id` y `name`) están activados, se considerarán duplicadas las secuencias que tengan el mismo `id` **y** el mismo `name`.
        - Si solo uno de los parámetros está activado:
          - `id=True, name=False`: Elimina duplicados por `id`, conservando la primera aparición.
          - `id=False, name=True`: Elimina duplicados por `name`.
        - Si ambos parámetros están desactivados (`id=False, name=False`), se mostrará un mensaje de advertencia y no se realizará ninguna acción.
        - La eliminación de duplicados afecta solo a la lista de secuencias. Se recomienda reindexar posteriormente para actualizar los diccionarios.
        """

        # Si no se especifica ningún criterio, se muestra una advertencia y se sale del método
        if not id and not name:
            logging.info("No se ha seleccionado ningún criterio de duplicados. No se realizaron cambios.")
            return

        # Eliminación de duplicados considerando ambos criterios: `id` y `name`
        if id and name:
            # Crea un diccionario con clave (id, name), lo que garantiza unicidad por ambas propiedades
            unique = {(seq.id, seq.name): seq for seq in self.lista_secuencias}
            self.lista_secuencias = list(unique.values())

        # Eliminación de duplicados solo por `id`
        elif id:
            unique = {seq.id: seq for seq in self.lista_secuencias}  # Última ocurrencia prevalece
            self.lista_secuencias = list(unique.values())

        # Eliminación de duplicados solo por `name`
        elif name:
            unique = {seq.name: seq for seq in self.lista_secuencias}  # Última ocurrencia prevalece
            self.lista_secuencias = list(unique.values())

    def guardar_records(self, lista_paths: List[str], lista_format: List[str], output_name: str):
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

    def eliminar_secuencias_by_id(self, id: str = None, all: bool = False) -> int:
        """
        Elimina secuencias de la lista y actualiza los índices.

        Parámetros: \n
        - id (str, opcional): Identificador de la secuencia a eliminar.
        - all (bool, opcional): Si es `True`, elimina todas las secuencias.

        Retorna: \n
        - int: Número de secuencias eliminadas.

        Excepciones: \n
        - ValueError: Si se proporciona un ID junto con `all=True`.
        - KeyError: Si el ID no se encuentra en la lista de secuencias.
        """

        # Validación de parámetros
        if id and all:
            raise ValueError("Si el parámetro 'all' es verdadero, no se debe introducir un ID.")

        # Eliminar todas las secuencias
        if all:
            num_eliminadas = len(self.lista_secuencias)
            self.lista_secuencias.clear()
            self.index_sequences_by_id_name(force=True)
            logging.info(f"Se eliminaron {num_eliminadas} secuencias.")
            return num_eliminadas

        # Buscar y eliminar una secuencia específica
        if id:
            record_a_eliminar = [record for record in self.lista_secuencias if record.id == id]

            # Si el ID no existe, lanzar excepción
            if not record_a_eliminar:
                raise KeyError(f"ID '{id}' no encontrado en la lista de secuencias.")

            # Nueva lista sin la secuencia eliminada (más eficiente que `remove()`)
            self.lista_secuencias = [record for record in self.lista_secuencias if record.id != id]

            # Reindexar después de la eliminación
            self.index_sequences_by_id_name(force=True)

            logging.info(f"Se eliminó la secuencia con ID '{id}'.")
            return len(record_a_eliminar)

        return 0  # Si no se especifica un ID ni `all`, no se realiza ninguna acción.

    def asignar_cds(self, seq_id: str, inicio: int, final: int):
        """
           Asigna manualmente una región codificante (CDS) a una secuencia.

           Parámetros:
           - seq_id (str): ID de la secuencia.
           - inicio (int): Posición de inicio del CDS.
           - final (int): Posición de fin del CDS.

           Retorna: \n
           None (En su lugar modifica el diccionario cds_dict)

           Excepciones: \n
           - ValueError: Valores de inicio o final invalidos o menores de 0.
           - KeyError: El ID no existe en la lista de secuencias.
           """

        if seq_id not in self.index_by_id.keys():
            raise KeyError(f"{seq_id} no existe en la lista de secuencias.")
        if inicio > final or inicio < 0 or final < 0:
            raise ValueError("Coordenadas de cds invalidas")

        self.cds_dict[seq_id] = (inicio, final)
        logging.info(f"CDS asignado a {seq_id}: ({inicio}, {final})")

    def obtener_secuencia_sin_exones(self, seq_id: str) -> str:
        """
        Obtiene la secuencia de ADN excluyendo los intrones, basada en los exones o el CDS si están disponibles.

        Parámetros:
        - seq_id (str): ID de la secuencia.

        Retorna:
        - str: Secuencia de ADN sin intrones (si hay exones) o la secuencia CDS si está definida.
               Si no hay información de exones ni CDS, devuelve la secuencia completa.

        Excepciones:
        - ValueError: Si el ID de la secuencia no existe en `index_by_id` o es inválido.

        Funcionamiento:
        1. Verifica si `seq_id` existe en `index_by_id`. Si no, lanza un error.
        2. Obtiene la secuencia completa usando `get_sequence_by_id(seq_id)`.
        3. Si hay información de exones en `exones_dict`, reconstruye la secuencia sin intrones.
        4. Si no hay exones pero hay un CDS en `cds_dict`, usa el CDS como referencia.
        5. Si no hay ninguna anotación, devuelve la secuencia completa sin cambios.
        """

        # Validación: Verificar si el ID existe en el diccionario de secuencias
        if not seq_id or seq_id not in self.index_by_id:
            raise ValueError(f"❌ ID de secuencia '{seq_id}' no existe en la lista de secuencias.")

        # Obtener la secuencia completa desde el diccionario de secuencias
        secuencia = self.get_sequence_by_id(seq_id)

        # Si la secuencia tiene exones registrados, unir solo las partes correspondientes
        if seq_id in self.exones_dict:
            seq_sin_exones =  "".join(str(secuencia.seq)[start:end + 1] for start, end in self.exones_dict[seq_id])
            secuencia.seq = Seq(seq_sin_exones)
            return secuencia.seq

        # Si la secuencia no tiene exones pero tiene una región CDS, usar solo el CDS
        if seq_id in self.cds_dict:
            start, end = self.cds_dict[seq_id]
            secuencia.seq = secuencia.seq[start:end + 1]
            return secuencia.seq

        # Si no hay información de exones ni CDS, devolver la secuencia completa
        return secuencia.seq

    def filtrar_secuencias_por_longitud(self, min_len: int, max_len: int, reindex: bool = True) -> int:
        """
        Filtra las secuencias de la lista interna `lista_secuencias` según un rango de longitud especificado.

        Parámetros: \n
        - min_len (int): Longitud mínima permitida para las secuencias.
        - max_len (int): Longitud máxima permitida para las secuencias.

        Retorna: \n
        - int: Número de secuencias que quedaron después de la filtración.

        Excepciones: \n
        - ValueError: Si los valores de `min_len` o `max_len` no son válidos o no cumplen con las condiciones requeridas.
        - TypeError: Si los argumentos proporcionados no son del tipo esperado.

        Funcionamiento: \n
        1. **Validación previa**: Se comprueba que `min_len` y `max_len` sean enteros positivos y que `min_len < max_len`.
        2. **Filtrado**: Se llama a `filter_sequences_by_length(min_len, max_len, self.lista_secuencias)`,
           que filtra las secuencias dentro del rango de longitud especificado.
        3. **Manejo de errores**: Captura posibles excepciones (`ValueError` y `TypeError`) y las imprime.
        4. **Retorno**: Devuelve el número de secuencias que cumplen con el filtro.

        Notas: \n
        - Si `min_len` es demasiado alto o `max_len` demasiado bajo, `self.lista_secuencias` podría quedar vacía.
        - Se recomienda reindexar después de la filtración si se requiere acceso rápido por `id` o `name`.
        - Este método **modifica `self.lista_secuencias` en su lugar**, por lo que las secuencias que no cumplan
          con el filtro se eliminan permanentemente.
        """

        # 🔍 Validación de entrada antes de procesar
        if not isinstance(min_len, int) or not isinstance(max_len,
                                                          int) or min_len < 0 or max_len < 0 or min_len > max_len:
            raise ValueError(
                "Los valores de min_len y max_len deben ser enteros positivos y cumplir min_len < max_len.")

        try:
            # Filtrar secuencias dentro del rango especificado
            self.lista_secuencias = list(filter_sequences_by_length(min_len, max_len, self.lista_secuencias))

        except ValueError as v:
            # Captura errores de valores inválidos
            logging.error(f"ValueError: {v}")

        except TypeError as t:
            # Captura errores de tipos de datos incorrectos
            logging.error(f"TypeError: {t}")

        if reindex:
            self.index_sequences_by_id_name(force=True)

        # Retornar el número de secuencias restantes después de la filtración
        return len(self.lista_secuencias)

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


    def convertir_formato(self, formato_salida: str, output_name: str):
        """
        Convierte las secuencias cargadas a un nuevo formato y las guarda en un archivo.

        Parámetros:
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
        if not self.lista_secuencias:
            raise ValueError("No hay secuencias cargadas para convertir.")

        # Asegurar que el nombre del archivo tenga la extensión correcta
        extension_dict = {"fasta": ".fasta", "genbank": ".gb", "fastq": ".fastq"}
        extension = extension_dict[formato_salida]
        output_name = output_name.strip()
        if not output_name.endswith(extension):
            output_name += extension

        # Solución para GenBank: Añadir molecule_type si falta
        if formato_salida == "genbank":
            for seq in self.lista_secuencias:
                if "molecule_type" not in seq.annotations:
                    if str(seq.seq).upper().find("U") == -1:
                        seq.annotations["molecule_type"] = "DNA"
                    else:
                        seq.annotations["molecule_type"] = "RNA"

        # Solución para FASTQ: Añadir valores de calidad ficticios si faltan
        if formato_salida == "fastq":
            for seq in self.lista_secuencias:
                if "phred_quality" not in seq.letter_annotations:
                    seq.letter_annotations["phred_quality"] = [40] * len(seq.seq)  # Calidad máxima ficticia

        # Guardar las secuencias en el formato especificado
        SeqIO.write(self.lista_secuencias, f"output/{output_name}", formato_salida)

        logging.info(f"Archivo guardado como '{output_name}' en formato '{formato_salida}'.")

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

    def traducir_secuencias(self, table: str or int = "Standard", to_stop: bool = False) -> Tuple[
        Dict[str, str], Dict[str, str]]:
        """
        Traduce las secuencias de ADN en `lista_secuencias` a proteínas usando una tabla de traducción específica.

        Parámetros:
        - table (str o int, opcional): Nombre o ID de la tabla de traducción genética.
          Por defecto, usa la tabla "Standard".
        - to_stop (bool, opcional): Si es `True`, la traducción se detiene en el primer codón de STOP.

        Retorna:
        - Tuple[Dict[str, str], Dict[str, str]]:
          1️. Diccionario con las secuencias traducidas `{ID_secuencia: Proteína}`.
          2️. Diccionario con secuencias que no se pudieron traducir `{ID_secuencia: Mensaje de error}`.

        Excepciones:
        - ValueError: Si `lista_secuencias` está vacía.
        - KeyError: Si la tabla de traducción no es válida.
        """

        # Verificar que haya secuencias
        if not self.lista_secuencias:
            raise ValueError("❌ Lista de secuencias vacía.")

        # Validar tabla de traducción
        try:
            if isinstance(table, str):
                _ = CodonTable.unambiguous_dna_by_name[table]  # Validamos existencia
            elif isinstance(table, int):
                _ = CodonTable.unambiguous_dna_by_id[table]
            else:
                raise TypeError("El parámetro table debe ser un nombre o un ID válido.")
        except KeyError:
            raise KeyError(
                f"Tabla de traducción '{table}' no encontrada. Tablas disponibles: {CodonTable.unambiguous_dna_by_id.keys()}")

        diccionario_traduccion = {}
        errores = {}

        # Traducir secuencias
        for seq_id, record in self.index_by_id.items():

            secuencia_madura = self.obtener_secuencia_sin_exones(seq_id)

            if len(secuencia_madura) % 3 != 0:
                errores[seq_id] = f"Longitud incorrecta ({len(secuencia_madura)} nucleótidos). No se puede traducir."
                continue  # Saltar a la siguiente secuencia

            try:
                traduccion = secuencia_madura.translate(table = table, to_stop = to_stop)
                diccionario_traduccion[seq_id] = traduccion
            except Exception as e:
                errores[seq_id] = f"Error en la traducción: {e}"

        # Mostrar errores encontrados (opcional)
        if errores:
            logging.info("\n Errores encontrados en las siguientes secuencias:")
            for id_seq, error in errores.items():
                logging.info(f" {id_seq}: {error}")

        return diccionario_traduccion, errores


def menu():
    """Muestra las opciones del menú al usuario."""
    print("\nMenú Principal")
    print("1️⃣  Cargar secuencias desde un archivo")
    print("2️⃣  Mostrar secuencias cargadas")
    print("3️⃣  Buscar secuencia por ID o nombre")
    print("4️⃣  Eliminar secuencias o duplicados")
    print("5️⃣  Filtrar secuencias por longitud")
    print("6️⃣  Buscar motivos en las secuencias")
    print("7️⃣  Mostrar y asignar CDS/Exones")
    print("8️⃣  Traducir secuencias")
    print("9️⃣  Obtener estadísticas de secuencias")
    print("🔟  Exportar secuencias en otro formato")
    print("📥  11 - Combinar múltiples archivos en un solo FASTA")
    print("0️⃣  Salir")


def cargar_archivo(pipeline):
    """Función para cargar un archivo de secuencias."""
    file_path = input("Ingresa la ruta del archivo: ")
    file_format = input("Ingresa el formato del archivo (fasta, genbank, fastq): ").lower()
    pipeline.load_sequences(file_path, file_format)
    print(f"Se cargaron {len(pipeline.lista_secuencias)} secuencias.")


def mostrar_secuencias(pipeline):
    """Muestra todas las secuencias cargadas."""
    if not pipeline.lista_secuencias:
        print("No hay secuencias cargadas.")
        return
    print(pipeline)


def buscar_secuencia(pipeline):
    """Busca una secuencia por ID o nombre."""
    criterio = input("Buscar por (id/nombre): ").lower()
    valor = input("Ingresa el ID o nombre: ")
    if criterio == "id":
        secuencia = pipeline.get_sequence_by_id(valor)
    elif criterio == "nombre":
        secuencia = pipeline.get_sequence_by_name(valor)
    else:
        print("Opción inválida.")
        return
    if secuencia:
        print(f"Secuencia encontrada:\n{secuencia}")
    else:
        print("Secuencia no encontrada.")


def eliminar_secuencias(pipeline):
    """Elimina secuencias por ID o todas."""
    opcion = input("Eliminar una secuencia (id) o todas (all): ").lower()
    if opcion == "id":
        seq_id = input("Ingresa el ID de la secuencia: ")
        eliminadas = pipeline.eliminar_secuencias_by_id(seq_id)
    elif opcion == "all":
        eliminadas = pipeline.eliminar_secuencias_by_id(all=True)
    else:
        print("Opción inválida.")
        return
    print(f"{eliminadas} secuencia(s) eliminada(s).")


def filtrar_secuencias(pipeline):
    """Filtra secuencias según su longitud."""
    try:
        min_len = int(input("Ingresa la longitud mínima: "))
        max_len = int(input("Ingresa la longitud máxima: "))
        cantidad = pipeline.filtrar_secuencias_por_longitud(min_len, max_len)
        print(f"Se filtraron {cantidad} secuencias en el rango especificado.")
    except ValueError:
        print("Ingresa valores numéricos válidos.")


def buscar_motivos(pipeline):
    """Busca motivos dentro de las secuencias."""
    motivos = input("Ingresa los motivos a buscar (separados por comas): ").split(",")
    resultados = pipeline.buscar_motivos(motivos)
    for seq_id, coincidencias in resultados.items():
        print(f"En la secuencia {seq_id}:")
        for motivo, posiciones in coincidencias.items():
            print(f"{motivo}: {posiciones}")


def mostrar_y_asignar_cds(pipeline):
    """Muestra CDS/exones y permite asignar manualmente un CDS."""
    seq_id = input("Ingresa el ID de la secuencia: ")
    if seq_id in pipeline.cds_dict:
        print(f"CDS actual: {pipeline.cds_dict[seq_id]}")
    if seq_id in pipeline.exones_dict:
        print(f"Exones actuales: {pipeline.exones_dict[seq_id]}")

    opcion = input("¿Quieres asignar un nuevo CDS? (y/n): ").lower()
    if opcion == "y":
        try:
            inicio = int(input("Ingresa la posición de inicio: "))
            final = int(input("Ingresa la posición final: "))
            pipeline.asignar_cds(seq_id, inicio, final)
            print(f"Nuevo CDS asignado: ({inicio}, {final})")
        except ValueError:
            print("Ingresa valores numéricos válidos.")


def traducir_secuencias(pipeline):
    """Traduce secuencias de ADN a proteínas."""
    traducciones, errores = pipeline.traducir_secuencias()
    for seq_id, proteina in traducciones.items():
        print(f"{seq_id}: {proteina}")
    for seq_id, error in errores.items():
        print(f"{seq_id}: {error}")


def obtener_estadisticas(pipeline):
    """Muestra estadísticas de las secuencias cargadas."""
    estadisticas = pipeline.obtener_estadisticas()
    for key, value in estadisticas.items():
        print(f"{key}: {value}")


def exportar_secuencias(pipeline):
    """Convierte y guarda secuencias en un nuevo formato."""
    formato = input("Formato de salida (fasta, genbank, fastq): ").lower()
    nombre = input("Nombre del archivo de salida: ")
    pipeline.convertir_formato(formato, nombre)
    print(f"Archivo guardado como {nombre}.{formato}")


def guardar_records_prueba(pipeline):
    """Carga secuencias desde múltiples archivos y guarda todas en un único archivo FASTA en 'data/'."""
    print("📂 Vas a combinar múltiples archivos de secuencias en un solo archivo FASTA.")

    # Pedir al usuario las rutas de los archivos
    lista_paths = input("📝 Ingresa las rutas de los archivos (separadas por comas): ").split(",")

    # Eliminar espacios en blanco alrededor de cada ruta
    lista_paths = [path.strip() for path in lista_paths]

    #Añadimos el directorio de donde provienven las secuencias
    lista_paths = ["data/" + path for path in lista_paths]

    # Verificar que los archivos existen
    archivos_no_encontrados = [path for path in lista_paths if not os.path.exists(path)]
    if archivos_no_encontrados:
        print(f"⚠️ Los siguientes archivos no existen: {', '.join(archivos_no_encontrados)}")
        return

    # Pedir formatos para cada archivo
    print("📄 Ahora ingresa los formatos de cada archivo en el mismo orden (fasta, genbank, fastq)")
    lista_format = input("➡️ Formatos (separados por comas): ").lower().split(",")

    # Eliminar espacios en blanco alrededor de cada formato
    lista_format = [fmt.strip() for fmt in lista_format]

    # Verificar que el número de archivos y formatos coincida
    if len(lista_paths) != len(lista_format):
        print("❌ Error: La cantidad de rutas y formatos no coincide. Intenta de nuevo.")
        return

    # Pedir nombre del archivo de salida
    output_name = input("📁 Ingresa el nombre del archivo de salida (sin extensión): ").strip()


    # Ejecutar guardar_records
    pipeline.guardar_records(lista_paths, lista_format, output_name)

    print(f"✅ Archivo combinado guardado en 'data/{output_name}.fasta'.")


def main():
    pipeline = Pipeline()
    while True:
        menu()
        opcion = input("Ingresa el número de la opción: ")
        if opcion == "1":
            cargar_archivo(pipeline)
        elif opcion == "2":
            mostrar_secuencias(pipeline)
        elif opcion == "3":
            buscar_secuencia(pipeline)
        elif opcion == "4":
            eliminar_secuencias(pipeline)
        elif opcion == "5":
            filtrar_secuencias(pipeline)
        elif opcion == "6":
            buscar_motivos(pipeline)
        elif opcion == "7":
            mostrar_y_asignar_cds(pipeline)
        elif opcion == "8":
            traducir_secuencias(pipeline)
        elif opcion == "9":
            obtener_estadisticas(pipeline)
        elif opcion == "10":
            exportar_secuencias(pipeline)
        elif opcion == "11":
            guardar_records_prueba(pipeline)
        elif opcion == "0":
            print("Saliendo del programa...")
            break
        else:
            print("Opción no válida, intenta nuevamente.")


if __name__ == "__main__":
    main()










