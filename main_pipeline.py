from pipeline import Pipeline

def ejecutar_pipeline():
    pipeline = Pipeline()
    pipeline.load_sequences("data/entrada.fasta")  # o desde argumentos

    # Opcional: retrotraducir
    pipeline.retrotraducir_a_adn(["prot1", "prot2"])

    # Ejecutar alineamiento
    pipeline.ejecutar_smith_waterman("blosum62", 2, "resultados.txt")

    # Otras operaciones...

if __name__ == "__main__":
    ejecutar_pipeline()
