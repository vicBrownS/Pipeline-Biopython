# 🔬 **Pipeline de Manejo de Secuencias Biológicas con Biopython**  

🚀 **Un pipeline en Python para la gestión, análisis y procesamiento de secuencias biológicas (FASTA, GenBank, FASTQ) utilizando Biopython.**  

---

## 📌 **Características principales**  
✅ Carga, almacenamiento y exportación de secuencias en distintos formatos (FASTA, GenBank, FASTQ).  
✅ Búsqueda de secuencias por ID o nombre, eliminación de duplicados y filtrado por longitud.  
✅ **Procesamiento avanzado de anotaciones**: extracción de regiones codificantes (CDS, exones e intrones).  
✅ **Traducción de secuencias de ADN a proteínas**, respetando marcos de lectura y tablas genéticas.  
✅ **Interfaz interactiva** con menú para facilitar su uso sin conocimientos previos en Biopython.  

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
*(El archivo `requirements.txt` debe contener Biopython y otras dependencias necesarias.)*

---

## 🚀 **Uso del pipeline**  
Ejecuta el `main.py` para acceder al menú interactivo:  
```bash
python main.py
```
### 📂 **Opciones disponibles en el menú**  
1️⃣ **Cargar secuencias** desde un archivo (FASTA, GenBank, FASTQ).  
2️⃣ **Mostrar secuencias cargadas** con ID, nombre y longitud.  
3️⃣ **Buscar secuencias** por ID o nombre.  
4️⃣ **Eliminar secuencias** específicas o eliminar duplicados.  
5️⃣ **Filtrar secuencias por longitud**.  
6️⃣ **Buscar motivos** dentro de las secuencias.  
7️⃣ **Mostrar y asignar CDS/exones** en secuencias GenBank.  
8️⃣ **Traducir secuencias** de ADN a proteínas.  
9️⃣ **Obtener estadísticas** como longitud media y contenido GC.  
🔟 **Exportar secuencias** a otro formato.  
📥 **Combinar múltiples archivos** en un solo FASTA.  

---

## 📄 **Ejemplo de uso**  
### 🧬 **Cargar un archivo y traducir secuencias**  
```bash
python main.py
```
*(En el menú, elige la opción 1 para cargar un archivo y la opción 8 para traducir las secuencias.)*

---

## 🛠 **Tecnologías utilizadas**  
- **Python**  
- **Biopython**  
- **Logging** para registro de eventos  
- **argparse** para línea de comandos  

---

## 📢 **Contribuciones**  
¡Las contribuciones son bienvenidas! Si quieres mejorar este proyecto, abre un **issue** o envía un **pull request**.  

📧 Contacto: https://www.linkedin.com/in/victor-brown-47050533a?utm_source=share&utm_campaign=share_via&utm_content=profile&utm_medium=android_app  

