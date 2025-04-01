import utils.traducciones
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from collections import Counter
from typing import List, Dict, Tuple
from utils.alineamiento import smith_waterman, medir_rendimiento_custom, medir_rendimiento_biopython, guardar_resultados_alineamientos
from utils.io_utils import guardar_records, convertir_formato, load_sequences, get_sequences_by_id
from utils.traducciones import traducir_secuencias, retrotraducir_a_adn
from utils.analisis import filter_sequences_by_length, obtener_estadisticas, buscar_motivos
import logging
import itertools

class Pipeline:
    """
    Clase que gestiona el procesamiento de secuencias biológicas (ADN, ARN y proteínas) mediante una arquitectura modular.

    **Funcionalidades principales:**

    1. **Carga y Gestión de Secuencias:**
       - `load_sequences(file_path, file_format)`: Carga secuencias desde archivos FASTA, GenBank o FASTQ.
       - `add_sequences(sequences, reindex=True)`: Añade nuevas secuencias manualmente.
       - `eliminar_secuencias_by_id(seq_id, all=False)`: Elimina secuencias específicas o todas.
       - `eliminar_duplicados(id=True, name=True)`: Elimina duplicados por ID o nombre.
       - `index_sequences_by_id_name(force=False)`: Indexa secuencias para búsqueda rápida.
       - `get_sequence_by_id(seq_id)` / `get_sequence_by_name(name)`: Accede a secuencias por ID o nombre.

    2. **Procesamiento de Secuencias:**
       - `filtrar_secuencias_por_longitud(...)`: Filtra secuencias por longitud mínima/máxima.
       - `buscar_motivos(motivos)`: Encuentra posiciones de motivos en las secuencias.
       - `asignar_cds(seq_id, inicio, final)`: Asigna regiones codificantes (CDS) manualmente.
       - `obtener_secuencia_sin_exones(seq_id)`: Extrae la secuencia sin exones definidos.

    3. **Traducción y Retrotraducción:**
       - `traducir_secuencias(table="Standard", to_stop=False)`: Traduce ADN a proteína usando la tabla genética especificada.
       - `retrotraducir_a_adn(ids)`: Convierte secuencias proteicas a secuencias de ADN mediante codones canónicos.

    4. **Alineamiento de Secuencias (Smith-Waterman):**
       - `ejecutar_smith_waterman(matrix_name, gap_penalty, guardar_en=None)`: Aplica el algoritmo de alineamiento local Smith-Waterman entre todas las combinaciones de secuencias.
       - Usa una matriz de sustitución (ej. BLOSUM62) y penalización por gap.
       - Permite guardar los resultados en un archivo `.txt`.
       - Comparación de rendimiento disponible con Biopython (`medir_rendimiento_*` desde `utils.alineamiento`).

    5. **Conversión y Guardado de Archivos:**
       - `convertir_formato(formato_salida, output_name)`: Convierte las secuencias actuales a otro formato (FASTA, GenBank, etc.).
       - `guardar_records(lista_paths, lista_format, output_name)`: Carga varios archivos de entrada y los fusiona en uno solo.

    6. **Análisis General:**
       - `obtener_estadisticas()`: Calcula estadísticas como longitud media, mínima, máxima, etc.

    **Notas:**
    - Toda la lógica pesada está modularizada en `utils/` (traducciones, alineamientos, análisis, E/S).
    - Las secuencias se gestionan mediante `SeqRecord` de Biopython.
    - Se puede utilizar el pipeline desde menú interactivo (`menu_pipeline.py`) o desde un flujo automatizado (`main_pipeline.py`).
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
        """Utiliza la función guardar_records del módulo io_utils"""
        guardar_records(lista_paths, lista_format, output_name)

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
        Utiliza la funcion de obtener_estadísticas del modulo analisis.py
        """
        return obtener_estadisticas(self)


    def convertir_formato(self, formato_salida: str, output_name: str):
        """
        Utiliza la funcion convertir formato del módulo io_utils
        """
        convertir_formato(self.lista_secuencias, formato_salida, output_name)

    def buscar_motivos(self, motivos: List[str]) -> Dict[str, Dict[str, List[int]]]:
        """
        Utiliza la función de buscar_motivos del modulo de analisis.py
        """
        diccionario_motivos = buscar_motivos(self,motivos)
        return diccionario_motivos

    def traducir_secuencias(self, table: str or int = "Standard", to_stop: bool = False) -> Tuple[
        Dict[str, str], Dict[str, str]]:
        """
        Utiliza la funcion de traducir secuencias del módulo traducciones.py
        """
        proteinas, errores = traducir_secuencias(table, to_stop)
        return proteinas, errores


    def retrotraducir_a_adn(self, ids: list[str]):
        """
        Utiliza la funcion de retrotraducir_a_adn del módulo traducción
        """
        retrotraducir_a_adn(self, ids)

    def ejecutar_smith_waterman(self,matrix_name: str, gap_penalty: int, guardar_en: str = None):
        """
        Ejecuta el algoritmo Smith-Waterman sobre todas las combinaciones de pares de secuencias cargadas.

        Parámetros:
            matrix_name (str): Nombre de la matriz de sustitución (ej. 'blosum62').
            gap_penalty (int): Penalización por inserción/eliminación.
            guardar_en (str): Nombre del archivo donde guardar los resultados (opcional).

        Retorna:
            list: Resultados como lista de tuplas (id1, id2, aln1, aln2, score)
        """
        resultados = []
        for s1, s2 in itertools.combinations(self.lista_secuencias, 2):
            aln1, aln2, score = smith_waterman(str(s1.seq), str(s2.seq), matrix_name, gap_penalty)
            resultados.append((s1.id, s2.id, aln1, aln2, score))

        if guardar_en:
            guardar_resultados_alineamientos(resultados, guardar_en)
            print(f"✅ Resultados guardados en '{guardar_en}'")

        return resultados
