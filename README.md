
# 🔬 **Pipeline de Manejo de Secuencias Biológicas con Biopython**

🚀 **Un pipeline en Python modularizado para la gestión, análisis, traducción y alineamiento de secuencias biológicas (FASTA, GenBank, FASTQ) utilizando Biopython.**

---

## 📌 **Características principales**

✅ Carga, almacenamiento y exportación de secuencias en distintos formatos (FASTA, GenBank, FASTQ).  
✅ Búsqueda de secuencias por ID o nombre, eliminación de duplicados y filtrado por longitud.  
✅ Extracción de anotaciones genómicas: CDS, exones e intrones.  
✅ **Traducción y retrotraducción** entre ADN y proteínas.  
✅ **Alineamiento local (Smith-Waterman)** manual y comparativo con Biopython.  
✅ **Modularización completa** en carpetas `utils/` para mejorar la organización y escalabilidad.  
✅ **Interfaz interactiva y flujo automatizado** separados (`menu_pipeline.py` y `main_pipeline.py`).  

---

## ⚙️ **Instalación y requisitos**

### **1️⃣ Clonar el repositorio**
```bash
git clone https://github.com/vicBrownS/pipeline-biopython.git
cd pipeline-biopython
```

### **2️⃣ Instalar dependencias**
Asegúrate de tener **Python 3.8+** instalado y ejecuta:
```bash
pip install -r requirements.txt
```

---

## 🚀 **Modos de ejecución del pipeline**

### 🧪 **Modo interactivo (menú)**
Ideal para explorar y probar funcionalidades paso a paso.
```bash
python menu_pipeline.py
```

### 🧬 **Modo automatizado (flujo real)**
Ejecuta un flujo completo definido en código.
```bash
python main_pipeline.py
```

---

## 📂 **Opciones del menú interactivo**

1️⃣ Cargar secuencias desde archivo  
2️⃣ Mostrar secuencias cargadas  
3️⃣ Buscar por ID o nombre  
4️⃣ Eliminar secuencias o duplicados  
5️⃣ Filtrar por longitud  
6️⃣ Buscar motivos en secuencias  
7️⃣ Asignar CDS y mostrar regiones codificantes  
8️⃣ Traducir ADN a proteína  
9️⃣ Obtener estadísticas (longitud media, GC...)  
🔟 Exportar secuencias a otros formatos  
🔁 Fusionar archivos en uno solo  
1️⃣2️⃣ Alineamiento local Smith-Waterman entre secuencias  
1️⃣3️⃣ Comparar rendimiento (propio vs Biopython)  
1️⃣4️⃣ Retrotraducir proteínas a ADN  

---

## 🛠️ **Estructura modular del proyecto**

```
proyecto_pipeline/
├── pipeline.py              # Clase principal del pipeline
├── menu_pipeline.py         # Interfaz interactiva con input del usuario
├── main_pipeline.py         # Flujo automatizado de ejecución
├── utils/
│   ├── io_utils.py          # Entrada/salida de archivos
│   ├── traducciones.py      # Traducción y retrotraducción ADN ↔ proteína
│   ├── alineamiento.py      # Smith-Waterman + comparación de rendimiento
│   └── analisis.py          # Filtros, estadísticas, CDS, motivos
```

---

## 📄 **Ejemplo de uso básico**

```bash
python menu_pipeline.py
```
(En el menú, elige la opción 1 para cargar secuencias y la opción 8 para traducirlas)

---

## 🧪 **Alineamiento Smith-Waterman (propio vs Biopython)**

Puedes ejecutar alineamientos locales con tu propia implementación y compararlos con Biopython:

```bash
python menu_pipeline.py
# Opción 12: ejecutar alineamiento
# Opción 13: comparar tiempos de ejecución
```

---

## 🛠 **Tecnologías utilizadas**

- **Python 3.8+**
- **Biopython**
- **Logging** para trazabilidad del flujo
- **argparse / input** para interacción desde terminal

---

## 🤝 **Contribuciones**

¡Las contribuciones son bienvenidas!  
Puedes abrir un **issue** o enviar un **pull request** con mejoras, ejemplos o nuevas funcionalidades.

📧 Contacto: [LinkedIn](https://www.linkedin.com/in/victor-brown-47050533a?utm_source=share&utm_campaign=share_via&utm_content=profile&utm_medium=android_app)
