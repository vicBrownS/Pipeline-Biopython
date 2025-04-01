"""
Módulo de traducción y retrotraducción de secuencias.

Agrupa funciones relacionadas con la conversión entre secuencias de ADN y proteínas, manteniendo separado el procesamiento biológico básico del resto del pipeline.

Funciones incluidas:
- `traducir_secuencias`: Traduce múltiples secuencias desde la clase Pipeline.
- `retrotraducir_a_adn`: Retrotraduce múltiples secuencias gestionadas desde Pipeline.

Este módulo permite mantener separadas las operaciones moleculares comunes y facilita su mantenimiento y mejora.
"""

from typing import Dict
from Bio.Data import CodonTable
import logging
from Bio.SeqIO import SeqRecord
from Bio.Seq import Seq

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
            traduccion = secuencia_madura.translate(table=table, to_stop=to_stop)
            diccionario_traduccion[seq_id] = traduccion
        except Exception as e:
            errores[seq_id] = f"Error en la traducción: {e}"

    # Mostrar errores encontrados (opcional)
    if errores:
        logging.info("\n Errores encontrados en las siguientes secuencias:")
        for id_seq, error in errores.items():
            logging.info(f" {id_seq}: {error}")

    return diccionario_traduccion, errores


def retrotraducir_a_adn(self, ids: list[str]):
    """
    Convierte secuencias de aminoácidos a ADN (retrotraducción) para los IDs indicados.

    Parámetros:
        ids (list[str]): Lista de IDs de secuencias de proteínas a convertir.

    Agrega nuevas secuencias de tipo ADN a la lista con un ID modificado.
    """
    codones_por_aminoacido = {
        'A': 'GCT', 'R': 'CGT', 'N': 'AAT', 'D': 'GAT', 'C': 'TGT',
        'Q': 'CAA', 'E': 'GAA', 'G': 'GGT', 'H': 'CAT', 'I': 'ATT',
        'L': 'CTT', 'K': 'AAA', 'M': 'ATG', 'F': 'TTT', 'P': 'CCT',
        'S': 'TCT', 'T': 'ACT', 'W': 'TGG', 'Y': 'TAT', 'V': 'GTT',
        '*': 'TAA', 'X': "AAA"  # Stop codón
    }

    for seq_id in ids:
        if seq_id not in self.index_by_id:
            print(f"❌ ID '{seq_id}' no encontrado.")
            continue

        record = self.index_by_id[seq_id]
        aa_seq = str(record.seq)

        try:
            adn = "".join(codones_por_aminoacido[aa] for aa in aa_seq)
        except KeyError as e:
            print(f"⚠️ Aminoácido no reconocido: {e} en secuencia {seq_id}")
            continue

        new_id = f"{seq_id}_retro"
        new_record = SeqRecord(Seq(adn), id=new_id, name=new_id, description="Retrotraducido desde proteína")
        self.add_sequences(new_record)

        print(f"✅ Secuencia retrotraducida añadida: {new_id}")
