## 🎯 **¿Qué es CIENBIO Sondas qPCR?**

CIENBIO es la **primera herramienta open-source especializada** en el diseño automatizado de sondas qPCR para análisis de polimorfismos de nucleótido único (SNPs). Combina algoritmos termodinámicos avanzados con un flujo de trabajo simplificado para generar sondas discriminatorias de alta precisión.

### 🚀 **Características Principales**

- ✅ **Diseño automático** de sondas SNP de 18-22 nucleótidos
- ✅ **Análisis termodinámico completo** con parámetros SantaLucia 1998
- ✅ **Cálculo de ΔG** a múltiples temperaturas (25°C y 37°C)
- ✅ **Evaluación de discriminación alélica** con ΔΔG cuantificado
- ✅ **Procesamiento masivo** desde archivos Excel
- ✅ **Análisis de estabilidad** termodinámica automático
- ✅ **Reportes estadísticos** detallados
- ✅ **Logging completo** para debugging científico

---

## 🤖 **Transparencia en el Desarrollo**

### **Metodología de Desarrollo Humano-IA**

Este proyecto fue desarrollado mediante **colaboración humano-IA**, estableciendo un nuevo estándar de transparencia en herramientas científicas:

**🧠 Aportaciones Humanas:**
- Supervisión científica y validación de algoritmos
- Verificación de parámetros termodinámicos (SantaLucia 1998)
- Definición de casos de uso clínicos
- Testing con datos reales de laboratorio
- Control de calidad científica

**🤖 Asistencia de IA (Claude 4 Sonnet - ChatGPT (OpenAI)):**
- Generación de código Python
- Estructuración orientada a objetos
- Implementación de manejo de errores
- Optimización de performance
- Documentación técnica

### **Validación Científica**
- ✅ Parámetros nearest-neighbor verificados contra literatura
- ✅ Cálculos de Tm validados con herramientas establecidas
- ✅ Casos de prueba con secuencias conocidas
- ✅ Revisión por científicos especializados

> *"Creemos que la transparencia en el desarrollo científico es fundamental para la confianza y reproducibilidad. La IA acelera el desarrollo, pero el expertise científico humano garantiza la calidad."*

---

## 📥 **Instalación**

### **Requisitos del Sistema**
- Python 3.7 o superior
- Sistema operativo: Windows, macOS, Linux

### **Dependencias**
```bash
pip install pandas openpyxl tkinter logging pathlib typing dataclasses traceback
```

### **Instalación Rápida**
```bash
# Clonar el repositorio
git clone https://github.com/cienbio2024/Sondas.git
cd Sondas

# Instalar dependencias
pip install -r requirements.txt

# Ejecutar el programa
python "Diseño de Sondas (IA)_refactored.py"
```

## 🚀 **Guía de Uso**

### **1. Preparar Archivo Excel**

Crea un archivo Excel con las siguientes columnas obligatorias:

| ID | Secuencia      | Coordenada SNP | Alelo_Ref | Alelo_Alt |
|----|----------------|----------------|-----------|-----------|
| SNP1 | ATCGATCGATCG.|      25        |    A      |     T     |
| SNP2 | GCTAGCTAGCTA.|      30        |    G      |     C     |

### **2. Ejecutar el Programa**

```bash
python "Diseño de Sondas (IA)_refactored.py"
```

1. Se abrirá un diálogo para seleccionar tu archivo Excel
2. El programa procesará automáticamente todas las secuencias
3. Se generarán los archivos de resultados

### **3. Archivos de Salida**

- **`archivo_evaluado.xlsx`**: Resultados completos con todas las propiedades
- **`archivo_resumen.txt`**: Resumen estadístico detallado
- **`cienbio_debug.log`**: Log técnico para debugging

### **Ejemplo de Uso Programático**

```python
from CalculadorTermodinamico import CalculadorTermodinamico
from GeneradorSondas import GeneradorSondas

# Inicializar calculadora
calculadora = CalculadorTermodinamico()

# Calcular propiedades de una secuencia
props = calculadora.calcular_propiedades_termodinamicas("ATGCATGCATGC")
print(f"Tm: {props['tm']}°C")
print(f"ΔG(25°C): {props['dg_25']} kcal/mol")
print(f"GC content: {props['gc_content']}%")

# Generar sondas para SNP
generador = GeneradorSondas(calculadora)
sondas = generador.generar_sondas_snp(
    secuencia="ATGCATGCATGCATGC",
    snp_pos=8,
    alelo_ref="G",
    alelo_alt="A",
    seq_id="test_snp"
)
```

---

## 📋 **Resultados Detallados**

### **Propiedades Calculadas**

Para cada sonda generada, el programa calcula:

- **Tm (°C)**: Temperatura de melting usando nearest-neighbor
- **ΔH (kcal/mol)**: Entalpía de formación
- **ΔG (kcal/mol)**: Energía libre de Gibbs a 25°C y 37°C
- **ΔΔG (kcal/mol)**: Diferencia de energía libre entre alelos
- **GC content (%)**: Porcentaje de bases G+C
- **Estabilidad**: Clasificación automática basada en ΔG

### **Interpretación de Resultados**

- **ΔG < 0**: Formación de dúplex favorecida (sonda estable)
- **ΔG > 0**: Formación de dúplex desfavorecida (sonda inestable)
- **|ΔΔG| > 2 kcal/mol**: Excelente discriminación alélica
- **GC 40-60%**: Rango óptimo para especificidad

---

## ⚙️ **Algoritmos y Métodos**

### **Base Científica**
- **Parámetros termodinámicos**: SantaLucia et al. (1998)
- **Método nearest-neighbor**: Para cálculos de Tm y ΔG
- **Correcciones de sal**: Ajustes para concentración iónica
- **Análisis de estabilidad**: Evaluación energética a múltiples temperaturas

### **Ecuaciones Principales**

```
Tm = (1000 * ΔH) / (ΔS + R * ln(C)) - 273.15 + 16.6 * log10([Na+])
ΔG = ΔH - T * ΔS
```

Donde:
- **ΔH**: Entalpía (kcal/mol)
- **ΔS**: Entropía (cal/mol·K)
- **C**: Concentración de oligonucleótido
- **T**: Temperatura (K)

---

## 🧪 **Validación y Testing**

### **Casos de Prueba Incluidos**
- Secuencias conocidas con propiedades publicadas
- SNPs reales de bases de datos como dbSNP
- Comparación con calculadoras establecidas (OligoCalc, Primer3)

### **Control de Calidad**
```bash
# Ejecutar tests (cuando estén disponibles)
python -m pytest tests/

# Validar con secuencia conocida
python validate_sequence.py ATGCATGCATGC
```

## 🐛 **Troubleshooting**

### **Problemas Comunes**

**Error: "Columnas faltantes"**
```
Solución: Verificar que el Excel tenga las columnas: ID, Secuencia, Coordenada SNP, Alelo_Ref, Alelo_Alt
```

**Error: "Secuencia inválida"**
```
Solución: Asegurar que las secuencias contengan solo nucleótidos A, T, G, C
```

**No se generan sondas**
```
Solución: Revisar el archivo cienbio_debug.log para detalles específicos
```

### **Logs y Debugging**

El programa genera logs detallados en `cienbio_debug.log`:
```bash
# Ver logs en tiempo real
tail -f cienbio_debug.log

# Buscar errores específicos
grep "ERROR" cienbio_debug.log
```

---

## 🤝 **Contribuciones**

### **Cómo Contribuir**

1. **Fork** el repositorio
2. Crea una **rama** para tu feature (`git checkout -b feature/nueva-funcionalidad`)
3. **Commit** tus cambios (`git commit -am 'Agregar nueva funcionalidad'`)
4. **Push** a la rama (`git push origin feature/nueva-funcionalidad`)
5. Abre un **Pull Request**

### **Tipos de Contribuciones Buscadas**

- 🐛 **Bug fixes** y mejoras de estabilidad
- 📚 **Documentación** y tutoriales
- 🧪 **Tests** unitarios y de integración
- ⚡ **Optimizaciones** de performance
- 🔬 **Validación científica** con nuevos datasets

### **Guidelines de Desarrollo**

- Mantener **transparencia** sobre contribuciones humanas vs. IA
- Incluir **tests** para nueva funcionalidad
- **Documentar** cambios en algoritmos científicos
- Seguir **PEP 8** para estilo de código Python

---

## 📄 **Licencia**

Este proyecto está licenciado bajo la **Licencia MIT** - ver el archivo [LICENSE](LICENSE) para detalles.

```
MIT License

Copyright (c) 2024 CIENBIO

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software...
```

---

## 📞 **Contacto y Soporte**

### **Equipo CIENBIO**
- 🌐 **Website**: [www.cienbio.cl](https://www.cienbio.cl)
- 📧 **Email**: contacto@cienbio.cl
- 🐙 **GitHub**: [@cienbio2024](https://github.com/cienbio2024)

### **Reportar Issues**
- 🐛 **Bugs**: [GitHub Issues](https://github.com/cienbio2024/Sondas/issues)
- 💡 **Feature Requests**: [GitHub Discussions](https://github.com/cienbio2024/Sondas/discussions)

### **Comunidad Científica**
- 🔬 **Colaboraciones**: Abiertos a colaboraciones académicas
- 📊 **Validaciones**: Compartir datos de validación
- 🎓 **Educación**: Uso en cursos de bioinformática

---

## 🏆 **Reconocimientos**

### **Algoritmos y Referencias Científicas**
- SantaLucia, J. (1998). *A unified view of polymer, dumbbell, and oligonucleotide DNA nearest-neighbor thermodynamics*. PNAS.
- Allawi, H.T. & SantaLucia, J. (1997). *Thermodynamics and NMR of internal G·T mismatches in DNA*. Biochemistry.

### **Transparencia en Desarrollo**
- **Claude 4 Sonnet** (Anthropic): Asistencia en generación de código
- **ChatGPT 4.0 (OpenAI): Asistencia en generación y revisión del código
- **Copilot Github: Revisión del código en VS Code
- **Supervisión científica humana**: Validación y control de calidad
- **Comunidad open-source**: Inspiración en transparencia y colaboración

### **Herramientas de Referencia**
- **Primer3**: Estándar académico para diseño de primers
- **OligoCalc**: Calculadora de propiedades termodinámicas
- **Biopython**: Biblioteca de bioinformática para Python

---

## 🔮 **Roadmap Futuro**

### **v2.2 - Próxima Versión**
- [ ] Interfaz web interactiva
- [ ] Soporte para análisis multiplex
- [ ] Integración con bases de datos SNP
- [ ] Visualizaciones gráficas de resultados

### **v2.3 - Características Avanzadas**
- [ ] Diseño de primers flanqueantes
- [ ] Análisis de especificidad genómica
- [ ] Exportación a formatos adicionales
- [ ] API REST para integración

### **v3.0 - Expansión de Funcionalidades**
- [ ] Soporte para InDels y CNVs
- [ ] Machine learning para optimización
- [ ] Interfaz gráfica completa
- [ ] Módulos de análisis estadístico avanzado

---
