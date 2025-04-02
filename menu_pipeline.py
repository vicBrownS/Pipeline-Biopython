from Pipeline import Pipeline
from utils.alineamiento import medir_rendimiento_custom, medir_rendimiento_biopython


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
    print("1️⃣2️⃣  Ejecutar Smith-Waterman sobre todas las secuencias")
    print("1️⃣3️⃣  Medir rendimiento de alineamientos (propio vs Biopython)")
    print("1️⃣4️⃣  Retrotraducir proteínas a ADN")

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
        elif opcion == "12":
            matrix = input("🔠 Matriz de sustitución (ej. BLOSUM62): ").strip()
            gap = int(input("➖ Penalización por gap (entero positivo): "))
            archivo = input("📁 Nombre del archivo de salida (opcional, pulsa Enter para omitir): ").strip()
            archivo = archivo if archivo else None
            pipeline.ejecutar_smith_waterman(matrix, gap, archivo)
        elif opcion == "13":
            matrix = input("🔠 Matriz de sustitución (ej. BLOSUM62): ").strip()
            gap = int(input("➖ Penalización por gap (entero positivo): "))

            print("\n⏱️  Medición con algoritmo propio:")
            medir_rendimiento_custom(pipeline.lista_secuencias, matrix, gap)

            print("\n🧬 Medición con PairwiseAligner de Biopython:")
            medir_rendimiento_biopython(pipeline.lista_secuencias, matrix)
        elif opcion == "14":
            ids = input("🧬 Ingresa los IDs de secuencias a retrotraducir (separados por comas): ").split(",")
            ids = [i.strip() for i in ids]
            pipeline.retrotraducir_a_adn(ids)
        elif opcion == "0":
            print("Saliendo del programa...")
            break
        else:
            print("Opción no válida, intenta nuevamente.")

if __name__ == "__main__":
    main()
