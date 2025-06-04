#!/usr/bin/env python3
"""
CIENBIO - Diseñador de Sondas qPCR v2.1
Desarrollado con asistencia de IA (Claude 4 Sonnet) bajo supervisión científica

DESARROLLO:
- Algoritmos termodinámicos: Supervisión científica humana
- Implementación Python: Asistencia de IA + revisión científica  
- Validación: Parámetros SantaLucia 1998 verificados manualmente
- Testing: Validación con casos de uso reales

TRANSPARENCIA:
Este software fue desarrollado mediante colaboración humano-IA:
- IA: Generación de código, estructuración, debugging
- Humano: Supervisión científica, validación de algoritmos, casos de uso

Autor científico: CIENBIO 
Asistencia técnica: Claude 4 Sonnet (Anthropic)
Licencia: MIT
Versión: 2.1
"""

import pandas as pd
import math
import tkinter as tk
from tkinter import filedialog, messagebox
import logging
from pathlib import Path
from typing import Dict, List, Tuple, Optional, Union
from dataclasses import dataclass
import traceback

# Configuración de logging más detallada
logging.basicConfig(
    level=logging.DEBUG,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.StreamHandler(),
        logging.FileHandler('cienbio_debug.log')
    ]
)
logger = logging.getLogger(__name__)


@dataclass
class SondaResult:
    """Clase para almacenar resultados de evaluación de sondas"""
    id_seq: str
    sonda_ref: str
    sonda_alt: str
    longitud: int
    inicio: int
    snp_central: int
    tm_ref: float
    tm_alt: float
    delta_tm: float
    dh_ref: float
    dh_alt: float
    dg_ref_25: float  # ΔG a 25°C
    dg_alt_25: float  # ΔG a 25°C
    dg_ref_37: float  # ΔG a 37°C
    dg_alt_37: float  # ΔG a 37°C
    gc_content: float  # Contenido GC %


class CalculadorTermodinamico:
    """Calculadora de propiedades termodinámicas para oligonucleótidos"""
    
    # Parámetros SantaLucia 1998 (ΔH, ΔS en kcal/mol y cal/mol·K)
    NN_PARAMS = {
        'AA': (-7.9, -22.2), 'TT': (-7.9, -22.2), 
        'AT': (-7.2, -20.4), 'TA': (-7.2, -21.3),
        'CA': (-8.5, -22.7), 'TG': (-8.5, -22.7), 
        'GT': (-8.4, -22.4), 'AC': (-8.4, -22.4),
        'CT': (-7.8, -21.0), 'AG': (-7.8, -21.0), 
        'GA': (-8.2, -22.2), 'TC': (-8.2, -22.2),
        'CG': (-10.6, -27.2), 'GC': (-9.8, -24.4), 
        'GG': (-8.0, -19.9), 'CC': (-8.0, -19.9)
    }
    
    def __init__(self, primer_conc: float = 500e-9, salt_conc: float = 50e-3):
        """
        Inicializa el calculador termodinámico
        
        Args:
            primer_conc: Concentración de primer en M (default: 500 nM)
            salt_conc: Concentración de sal en M (default: 50 mM)
        """
        self.primer_conc = primer_conc
        self.salt_conc = salt_conc
        self.R = 1.987  # Constante de gases en cal/mol·K
        logger.debug(f"Calculadora inicializada: primer_conc={primer_conc}, salt_conc={salt_conc}")
    
    def calcular_propiedades_termodinamicas(self, secuencia: str) -> Dict[str, float]:
        """
        Calcula todas las propiedades termodinámicas de la secuencia
        
        Args:
            secuencia: Secuencia de oligonucleótido
            
        Returns:
            Diccionario con Tm, ΔH, ΔG a diferentes temperaturas y %GC
        """
        try:
            secuencia_original = secuencia
            secuencia = self._limpiar_secuencia(secuencia)
            logger.debug(f"Secuencia original: '{secuencia_original}' -> Limpia: '{secuencia}'")
            
            if not self._validar_secuencia(secuencia):
                logger.error(f"Secuencia inválida: '{secuencia}' (original: '{secuencia_original}')")
                return {
                    'tm': float('nan'),
                    'dh': float('nan'),
                    'dg_25': float('nan'),
                    'dg_37': float('nan'),
                    'gc_content': float('nan')
                }
            
            dh_total, ds_total = self._calcular_parametros_nn(secuencia)
            logger.debug(f"Parámetros NN: ΔH={dh_total}, ΔS={ds_total}")
            
            # Correcciones de iniciación
            ds_total += -1.4  # Corrección entropía de iniciación
            dh_total += 0.2   # Corrección entalpía de iniciación
            
            # Cálculo de Tm usando ecuación SantaLucia
            tm = self._calcular_tm(dh_total, ds_total)
            
            # Cálculo de ΔG a diferentes temperaturas
            dg_25 = self._calcular_dg(dh_total, ds_total, 25.0)  # 25°C (estándar)
            dg_37 = self._calcular_dg(dh_total, ds_total, 37.0)  # 37°C (fisiológica)
            
            # Cálculo de contenido GC
            gc_content = self._calcular_gc_content(secuencia)
            
            logger.debug(f"Propiedades calculadas - Tm: {tm}°C, ΔH: {dh_total}, ΔG(25°C): {dg_25}, ΔG(37°C): {dg_37}, GC: {gc_content}%")
            
            return {
                'tm': round(tm, 2),
                'dh': round(dh_total, 2),
                'dg_25': round(dg_25, 2),
                'dg_37': round(dg_37, 2),
                'gc_content': round(gc_content, 2)
            }
            
        except Exception as e:
            logger.error(f"Error calculando propiedades para '{secuencia}': {e}")
            return {
                'tm': float('nan'),
                'dh': float('nan'),
                'dg_25': float('nan'),
                'dg_37': float('nan'),
                'gc_content': float('nan')
            }
    
    def _limpiar_secuencia(self, secuencia: str) -> str:
        """Limpia y normaliza la secuencia"""
        if pd.isna(secuencia):
            return ""
        return str(secuencia).upper().replace(" ", "").replace("\n", "").replace("\t", "").replace("\r", "")
    
    def _validar_secuencia(self, secuencia: str) -> bool:
        """Valida que la secuencia contenga solo nucleótidos válidos"""
        if not secuencia or len(secuencia) < 2:
            logger.debug(f"Secuencia muy corta: '{secuencia}' (longitud: {len(secuencia)})")
            return False
        
        bases_validas = set('ATGC')
        bases_secuencia = set(secuencia)
        bases_invalidas = bases_secuencia - bases_validas
        
        if bases_invalidas:
            logger.debug(f"Bases inválidas encontradas: {bases_invalidas} en secuencia: '{secuencia}'")
            return False
            
        return True
    
    def _calcular_parametros_nn(self, secuencia: str) -> Tuple[float, float]:
        """Calcula parámetros termodinámicos usando nearest neighbor"""
        dh_total = 0.0
        ds_total = 0.0
        
        for i in range(len(secuencia) - 1):
            dinucleotido = secuencia[i:i+2]
            delta_h, delta_s = self.NN_PARAMS.get(dinucleotido, (-7.0, -20.0))
            dh_total += delta_h
            ds_total += delta_s
            logger.debug(f"Dinucleótido {dinucleotido}: ΔH={delta_h}, ΔS={delta_s}")
        
        return dh_total, ds_total
    
    def _calcular_tm(self, dh: float, ds: float) -> float:
        """Calcula temperatura de melting"""
        try:
            # Verificar valores válidos para evitar errores matemáticos
            if ds == 0 or self.primer_conc <= 0 or self.salt_conc <= 0:
                raise ValueError(f"Valores inválidos: ΔS={ds}, primer_conc={self.primer_conc}, salt_conc={self.salt_conc}")
            
            tm = (1000 * dh) / (ds + (self.R * math.log(self.primer_conc))) - 273.15 + \
                 16.6 * math.log10(self.salt_conc)
            return tm
        except Exception as e:
            logger.error(f"Error calculando Tm: ΔH={dh}, ΔS={ds}, error={e}")
            raise
    
    def _calcular_dg(self, dh: float, ds: float, temperatura_celsius: float) -> float:
        """
        Calcula energía libre de Gibbs (ΔG) a una temperatura específica
        
        Args:
            dh: Entalpía en kcal/mol
            ds: Entropía en cal/mol·K
            temperatura_celsius: Temperatura en °C
            
        Returns:
            ΔG en kcal/mol
        """
        try:
            # Convertir temperatura a Kelvin
            temp_kelvin = temperatura_celsius + 273.15
            
            # ΔG = ΔH - T×ΔS
            # ΔH está en kcal/mol, ΔS en cal/mol·K, necesitamos convertir ΔS a kcal/mol·K
            ds_kcal = ds / 1000.0  # Convertir cal/mol·K a kcal/mol·K
            
            dg = dh - (temp_kelvin * ds_kcal)
            return dg
            
        except Exception as e:
            logger.error(f"Error calculando ΔG: ΔH={dh}, ΔS={ds}, T={temperatura_celsius}°C, error={e}")
            raise
    
    def _calcular_gc_content(self, secuencia: str) -> float:
        """
        Calcula el contenido de GC de la secuencia
        
        Args:
            secuencia: Secuencia de nucleótidos
            
        Returns:
            Porcentaje de contenido GC
        """
        if not secuencia:
            return 0.0
        
        gc_count = secuencia.count('G') + secuencia.count('C')
        total_bases = len(secuencia)
        
        return (gc_count / total_bases) * 100.0
    
    def calcular_tm_y_dh(self, secuencia: str) -> Tuple[float, float]:
        """
        Método de compatibilidad - calcula solo Tm y ΔH
        
        Args:
            secuencia: Secuencia de oligonucleótido
            
        Returns:
            Tupla con (Tm en °C, ΔH en kcal/mol)
        """
        props = self.calcular_propiedades_termodinamicas(secuencia)
        return props['tm'], props['dh']


class GeneradorSondas:
    """Generador de sondas para análisis de SNPs"""
    
    def __init__(self, calculadora: CalculadorTermodinamico):
        self.calculadora = calculadora
        logger.debug("GeneradorSondas inicializado")
    
    def generar_sondas_snp(self, secuencia: str, snp_pos: Union[int, str], alelo_ref: str, 
                          alelo_alt: str, seq_id: str, 
                          longitudes: Optional[List[int]] = None) -> List[SondaResult]:
        """
        Genera sondas para análisis de SNP
        
        Args:
            secuencia: Secuencia completa
            snp_pos: Posición del SNP (1-indexed)
            alelo_ref: Alelo de referencia
            alelo_alt: Alelo alternativo
            seq_id: Identificador de la secuencia
            longitudes: Lista de longitudes a evaluar
            
        Returns:
            Lista de resultados de sondas
        """
        logger.debug(f"Generando sondas para ID: {seq_id}")
        
        try:
            if longitudes is None:
                longitudes = list(range(18, 23))
            
            # Limpiar y validar datos de entrada
            secuencia = self.calculadora._limpiar_secuencia(secuencia)
            if not secuencia:
                logger.warning(f"Secuencia vacía para ID: {seq_id}")
                return []
            
            # Convertir posición SNP a entero
            try:
                snp_pos_int = int(float(str(snp_pos)))  # Maneja casos como "25.0"
            except (ValueError, TypeError):
                logger.error(f"Posición SNP inválida para ID {seq_id}: '{snp_pos}'")
                return []
            
            if snp_pos_int < 1 or snp_pos_int > len(secuencia):
                logger.error(f"Posición SNP fuera de rango para ID {seq_id}: {snp_pos_int} (secuencia longitud: {len(secuencia)})")
                return []
            
            snp_pos_0 = snp_pos_int - 1  # Convertir a 0-indexed
            alelo_ref = str(alelo_ref).upper().strip()
            alelo_alt = str(alelo_alt).upper().strip()
            
            logger.debug(f"Datos procesados - Secuencia: {len(secuencia)} bp, SNP pos: {snp_pos_int}, Alelos: {alelo_ref}->{alelo_alt}")
            
            # Verificar que el alelo de referencia coincida con la secuencia
            if secuencia[snp_pos_0] != alelo_ref:
                logger.warning(f"Alelo de referencia no coincide para ID {seq_id}: esperado '{alelo_ref}', encontrado '{secuencia[snp_pos_0]}'")
                # Continuar anyway, podría ser útil
            
            resultados = []
            
            for longitud in longitudes:
                logger.debug(f"Evaluando longitud {longitud} para ID {seq_id}")
                sondas = self._generar_sondas_longitud(
                    secuencia, snp_pos_0, alelo_ref, alelo_alt, longitud
                )
                
                for sonda_ref, sonda_alt, inicio in sondas:
                    try:
                        logger.debug(f"Calculando propiedades para sondas: ref='{sonda_ref}', alt='{sonda_alt}'")
                        
                        # Calcular propiedades para sonda de referencia
                        props_ref = self.calculadora.calcular_propiedades_termodinamicas(sonda_ref)
                        props_alt = self.calculadora.calcular_propiedades_termodinamicas(sonda_alt)
                        
                        resultado = SondaResult(
                            id_seq=seq_id,
                            sonda_ref=sonda_ref,
                            sonda_alt=sonda_alt,
                            longitud=longitud,
                            inicio=inicio + 1,  # Convertir a 1-indexed
                            snp_central=snp_pos_int,
                            tm_ref=props_ref['tm'],
                            tm_alt=props_alt['tm'],
                            delta_tm=round(props_alt['tm'] - props_ref['tm'], 2),
                            dh_ref=props_ref['dh'],
                            dh_alt=props_alt['dh'],
                            dg_ref_25=props_ref['dg_25'],
                            dg_alt_25=props_alt['dg_25'],
                            dg_ref_37=props_ref['dg_37'],
                            dg_alt_37=props_alt['dg_37'],
                            gc_content=props_ref['gc_content']  # GC content es igual para ambas sondas (solo cambia 1 nucleótido)
                        )
                        resultados.append(resultado)
                        logger.debug(f"Sonda generada exitosamente: Tm_ref={props_ref['tm']}, Tm_alt={props_alt['tm']}, ΔG_ref(25°C)={props_ref['dg_25']}")
                        
                    except Exception as e:
                        logger.error(f"Error calculando sonda {sonda_ref} para ID {seq_id}: {e}")
                        continue
            
            logger.info(f"Generadas {len(resultados)} sondas para ID: {seq_id}")
            return resultados
            
        except Exception as e:
            logger.error(f"Error general generando sondas para ID {seq_id}: {e}")
            logger.error(traceback.format_exc())
            return []
    
    def _generar_sondas_longitud(self, secuencia: str, snp_pos: int, 
                                alelo_ref: str, alelo_alt: str, 
                                longitud: int) -> List[Tuple[str, str, int]]:
        """Genera sondas de una longitud específica"""
        sondas = []
        half = longitud // 2
        inicio = max(snp_pos - half, 0)
        fin = min(inicio + longitud, len(secuencia))
        
        logger.debug(f"Generando sonda longitud {longitud}: inicio={inicio}, fin={fin}, SNP={snp_pos}")
        
        if fin - inicio == longitud and snp_pos in range(inicio, fin):
            sonda_ref = secuencia[inicio:fin]
            pos_relativa = snp_pos - inicio
            
            logger.debug(f"Sonda base: '{sonda_ref}', posición relativa SNP: {pos_relativa}")
            
            # Verificar que el alelo de referencia coincida (más permisivo)
            if pos_relativa < len(sonda_ref) and sonda_ref[pos_relativa] == alelo_ref:
                sonda_alt = (sonda_ref[:pos_relativa] + 
                           alelo_alt + 
                           sonda_ref[pos_relativa + 1:])
                sondas.append((sonda_ref, sonda_alt, inicio))
                logger.debug(f"Sonda válida generada: ref='{sonda_ref}', alt='{sonda_alt}'")
            else:
                logger.debug(f"Alelo no coincide en posición {pos_relativa}: esperado '{alelo_ref}', encontrado '{sonda_ref[pos_relativa] if pos_relativa < len(sonda_ref) else 'N/A'}'")
        else:
            logger.debug(f"Sonda no válida: longitud calculada={fin-inicio}, requerida={longitud}, SNP en rango={snp_pos in range(inicio, fin)}")
        
        return sondas


class ProcesadorExcel:
    """Procesador de archivos Excel con datos de SNPs"""
    
    COLUMNAS_REQUERIDAS = ["ID", "Secuencia", "Coordenada SNP", "Alelo_Ref", "Alelo_Alt"]
    
    def __init__(self, generador_sondas: GeneradorSondas):
        self.generador_sondas = generador_sondas
        logger.debug("ProcesadorExcel inicializado")
    
    def procesar_archivo(self, excel_path: str) -> pd.DataFrame:
        """
        Procesa archivo Excel y genera evaluación de sondas
        
        Args:
            excel_path: Ruta al archivo Excel
            
        Returns:
            DataFrame con resultados
        """
        try:
            logger.info(f"Leyendo archivo Excel: {excel_path}")
            df = pd.read_excel(excel_path)
            logger.info(f"Archivo leído exitosamente. Filas: {len(df)}, Columnas: {list(df.columns)}")
            
            self._validar_columnas(df)
            
            # Mostrar muestra de datos para debugging
            logger.debug("Primeras 3 filas del archivo:")
            for idx, fila in df.head(3).iterrows():
                logger.debug(f"Fila {idx}: {dict(fila)}")
            
            todos_resultados = []
            filas_procesadas = 0
            filas_con_error = 0
            
            for idx, fila in df.iterrows():
                try:
                    logger.debug(f"Procesando fila {idx + 1}/{len(df)}")
                    resultados_fila = self._procesar_fila(fila)
                    todos_resultados.extend(resultados_fila)
                    filas_procesadas += 1
                    
                    if resultados_fila:
                        logger.debug(f"Fila {idx + 1}: {len(resultados_fila)} sondas generadas")
                    else:
                        logger.warning(f"Fila {idx + 1}: No se generaron sondas")
                        
                except Exception as e:
                    filas_con_error += 1
                    logger.error(f"Error procesando fila {idx + 1}: {e}")
                    logger.error(f"Datos de la fila: {dict(fila)}")
                    continue
            
            logger.info(f"Procesamiento completado: {filas_procesadas} filas procesadas, {filas_con_error} errores, {len(todos_resultados)} sondas totales")
            
            if not todos_resultados:
                logger.warning("No se generaron resultados. Verificar datos de entrada.")
                return pd.DataFrame()
            
            return self._convertir_a_dataframe(todos_resultados)
            
        except Exception as e:
            logger.error(f"Error procesando archivo {excel_path}: {e}")
            logger.error(traceback.format_exc())
            raise
    
    def _validar_columnas(self, df: pd.DataFrame):
        """Valida que el DataFrame tenga las columnas requeridas"""
        columnas_faltantes = set(self.COLUMNAS_REQUERIDAS) - set(df.columns)
        if columnas_faltantes:
            error_msg = f"Columnas faltantes: {columnas_faltantes}. Columnas disponibles: {list(df.columns)}"
            logger.error(error_msg)
            raise ValueError(error_msg)
        logger.debug("Validación de columnas exitosa")
    
    def _procesar_fila(self, fila) -> List[SondaResult]:
        """Procesa una fila individual del DataFrame"""
        try:
            # Extraer datos con manejo de valores NaN
            seq_id = str(fila["ID"]) if not pd.isna(fila["ID"]) else f"seq_{hash(str(fila))}"
            secuencia = str(fila["Secuencia"]) if not pd.isna(fila["Secuencia"]) else ""
            snp_pos = fila["Coordenada SNP"]
            alelo_ref = str(fila["Alelo_Ref"]) if not pd.isna(fila["Alelo_Ref"]) else ""
            alelo_alt = str(fila["Alelo_Alt"]) if not pd.isna(fila["Alelo_Alt"]) else ""
            
            logger.debug(f"Procesando: ID={seq_id}, Secuencia={len(secuencia)} bp, SNP={snp_pos}, Alelos={alelo_ref}->{alelo_alt}")
            
            return self.generador_sondas.generar_sondas_snp(
                secuencia=secuencia,
                snp_pos=snp_pos,
                alelo_ref=alelo_ref,
                alelo_alt=alelo_alt,
                seq_id=seq_id
            )
        except Exception as e:
            logger.error(f"Error en _procesar_fila: {e}")
            raise
    
    def _convertir_a_dataframe(self, resultados: List[SondaResult]) -> pd.DataFrame:
        """Convierte lista de resultados a DataFrame"""
        if not resultados:
            logger.warning("Lista de resultados vacía")
            return pd.DataFrame()
        
        logger.info(f"Convirtiendo {len(resultados)} resultados a DataFrame")
        
        df_resultado = pd.DataFrame([
            {
                "ID": r.id_seq,
                "Sonda_ref": r.sonda_ref,
                "Sonda_alt": r.sonda_alt,
                "Longitud": r.longitud,
                "Inicio": r.inicio,
                "SNP Central": r.snp_central,
                "Tm_ref (°C)": r.tm_ref,
                "Tm_alt (°C)": r.tm_alt,
                "ΔTm (°C)": r.delta_tm,
                "ΔH_ref (kcal/mol)": r.dh_ref,
                "ΔH_alt (kcal/mol)": r.dh_alt,
                "ΔG_ref_25°C (kcal/mol)": r.dg_ref_25,
                "ΔG_alt_25°C (kcal/mol)": r.dg_alt_25,
                "ΔG_ref_37°C (kcal/mol)": r.dg_ref_37,
                "ΔG_alt_37°C (kcal/mol)": r.dg_alt_37,
                "ΔΔG_25°C (kcal/mol)": round(r.dg_alt_25 - r.dg_ref_25, 2),
                "ΔΔG_37°C (kcal/mol)": round(r.dg_alt_37 - r.dg_ref_37, 2),
                "GC_content (%)": r.gc_content,
                "Estabilidad_25°C": "Estable" if r.dg_ref_25 < 0 else "Inestable",
                "Estabilidad_37°C": "Estable" if r.dg_ref_37 < 0 else "Inestable"
            }
            for r in resultados
        ])
        
        logger.info(f"DataFrame creado exitosamente: {len(df_resultado)} filas")
        return df_resultado


class InterfazUsuario:
    """Interfaz gráfica para selección de archivos"""
    
    @staticmethod
    def seleccionar_archivo() -> Optional[str]:
        """
        Abre diálogo para seleccionar archivo Excel
        
        Returns:
            Ruta del archivo seleccionado o None
        """
        try:
            root = tk.Tk()
            root.withdraw()
            
            archivo = filedialog.askopenfilename(
                title="Selecciona archivo Excel con datos de SNPs",
                filetypes=[("Excel files", "*.xlsx"), ("Excel files", "*.xls"), ("All files", "*.*")]
            )
            
            root.destroy()
            return archivo if archivo else None
        except Exception as e:
            logger.error(f"Error en selección de archivo: {e}")
            return None
    
    @staticmethod
    def mostrar_mensaje(titulo: str, mensaje: str, tipo: str = "info"):
        """Muestra mensaje al usuario"""
        try:
            root = tk.Tk()
            root.withdraw()
            
            if tipo == "error":
                messagebox.showerror(titulo, mensaje)
            elif tipo == "warning":
                messagebox.showwarning(titulo, mensaje)
            else:
                messagebox.showinfo(titulo, mensaje)
            
            root.destroy()
        except Exception as e:
            logger.error(f"Error mostrando mensaje: {e}")
            print(f"{titulo}: {mensaje}")


def generar_resumen_estadistico(df: pd.DataFrame) -> str:
    """Genera un resumen estadístico de los resultados"""
    if df.empty:
        return "No hay datos para generar resumen."
    
    resumen = []
    resumen.append("📊 RESUMEN ESTADÍSTICO")
    resumen.append("=" * 40)
    resumen.append(f"Total de sondas generadas: {len(df)}")
    
    # Estadísticas de Tm
    resumen.append(f"\n🌡️ TEMPERATURA DE MELTING:")
    resumen.append(f"  Tm promedio (ref): {df['Tm_ref (°C)'].mean():.1f} ± {df['Tm_ref (°C)'].std():.1f} °C")
    resumen.append(f"  Tm promedio (alt): {df['Tm_alt (°C)'].mean():.1f} ± {df['Tm_alt (°C)'].std():.1f} °C")
    resumen.append(f"  ΔTm promedio: {df['ΔTm (°C)'].mean():.1f} ± {df['ΔTm (°C)'].std():.1f} °C")
    
    # Estadísticas de ΔG
    resumen.append(f"\n⚡ ENERGÍA LIBRE (ΔG):")
    resumen.append(f"  ΔG promedio a 25°C: {df['ΔG_ref_25°C (kcal/mol)'].mean():.2f} ± {df['ΔG_ref_25°C (kcal/mol)'].std():.2f} kcal/mol")
    resumen.append(f"  ΔG promedio a 37°C: {df['ΔG_ref_37°C (kcal/mol)'].mean():.2f} ± {df['ΔG_ref_37°C (kcal/mol)'].std():.2f} kcal/mol")
    resumen.append(f"  ΔΔG promedio a 25°C: {df['ΔΔG_25°C (kcal/mol)'].mean():.2f} ± {df['ΔΔG_25°C (kcal/mol)'].std():.2f} kcal/mol")
    resumen.append(f"  ΔΔG promedio a 37°C: {df['ΔΔG_37°C (kcal/mol)'].mean():.2f} ± {df['ΔΔG_37°C (kcal/mol)'].std():.2f} kcal/mol")
    
    # Estabilidad
    estables_25 = len(df[df['Estabilidad_25°C'] == 'Estable'])
    estables_37 = len(df[df['Estabilidad_37°C'] == 'Estable'])
    resumen.append(f"\n🔬 ESTABILIDAD:")
    resumen.append(f"  Sondas estables a 25°C: {estables_25}/{len(df)} ({estables_25/len(df)*100:.1f}%)")
    resumen.append(f"  Sondas estables a 37°C: {estables_37}/{len(df)} ({estables_37/len(df)*100:.1f}%)")
    
    # Contenido GC
    resumen.append(f"\n🧬 CONTENIDO GC:")
    resumen.append(f"  GC promedio: {df['GC_content (%)'].mean():.1f} ± {df['GC_content (%)'].std():.1f}%")
    resumen.append(f"  Rango GC: {df['GC_content (%)'].min():.1f}% - {df['GC_content (%)'].max():.1f}%")
    
    # Longitudes
    resumen.append(f"\n📏 LONGITUDES:")
    longitudes = df['Longitud'].value_counts().sort_index()
    for longitud, cantidad in longitudes.items():
        resumen.append(f"  {longitud} nt: {cantidad} sondas")
    
    return "\n".join(resumen)


def main():
    """Función principal del programa"""
    print("🧬 CIENBIO - Diseñador de Sondas qPCR v2.1")
    print("=" * 50)
    
    try:
        # Seleccionar archivo
        archivo_excel = InterfazUsuario.seleccionar_archivo()
        
        if not archivo_excel:
            print("❌ No se seleccionó ningún archivo.")
            return
        
        # Verificar que el archivo existe
        if not Path(archivo_excel).exists():
            error_msg = f"El archivo no existe: {archivo_excel}"
            logger.error(error_msg)
            InterfazUsuario.mostrar_mensaje("Error", error_msg, "error")
            return
        
        # Inicializar componentes
        logger.info("Inicializando componentes...")
        calculadora = CalculadorTermodinamico()
        generador = GeneradorSondas(calculadora)
        procesador = ProcesadorExcel(generador)
        
        print(f"📁 Procesando archivo: {Path(archivo_excel).name}")
        logger.info(f"Iniciando procesamiento de: {archivo_excel}")
        
        # Procesar archivo
        df_resultado = procesador.procesar_archivo(archivo_excel)
        
        if df_resultado.empty:
            mensaje_warning = "No se generaron resultados válidos. Revisa el archivo de log para más detalles."
            logger.warning(mensaje_warning)
            print(f"⚠️ {mensaje_warning}")
            InterfazUsuario.mostrar_mensaje("Advertencia", mensaje_warning, "warning")
            return
        
        # Guardar resultados
        archivo_salida = str(Path(archivo_excel).with_suffix('')) + "_evaluado.xlsx"
        logger.info(f"Guardando resultados en: {archivo_salida}")
        
        df_resultado.to_excel(archivo_salida, index=False)
        
        # Generar resumen estadístico
        resumen = generar_resumen_estadistico(df_resultado)
        logger.info("Resumen estadístico generado")
        
        # Guardar resumen en archivo de texto
        archivo_resumen = str(Path(archivo_excel).with_suffix('')) + "_resumen.txt"
        with open(archivo_resumen, 'w', encoding='utf-8') as f:
            f.write(resumen)
        
        # Mensaje de éxito
        mensaje_exito = (
            f"✅ Evaluación completada exitosamente!\n\n"
            f"📊 Sondas generadas: {len(df_resultado)}\n"
            f"💾 Archivo Excel: {Path(archivo_salida).name}\n"
            f"📄 Resumen estadístico: {Path(archivo_resumen).name}\n"
            f"📝 Log detallado: cienbio_debug.log\n\n"
            f"🔬 Propiedades calculadas:\n"
            f"  • Temperatura de melting (Tm)\n"
            f"  • Entalpía (ΔH)\n"
            f"  • Energía libre de Gibbs (ΔG) a 25°C y 37°C\n"
            f"  • Contenido GC (%)\n"
            f"  • Análisis de estabilidad termodinámica"
        )
        
        print(mensaje_exito)
        print("\n" + resumen)
        logger.info("Procesamiento completado exitosamente")
        InterfazUsuario.mostrar_mensaje("Éxito", 
            f"✅ Evaluación completada!\n\n"
            f"📊 {len(df_resultado)} sondas generadas\n"
            f"💾 Resultados: {Path(archivo_salida).name}\n"
            f"📄 Resumen: {Path(archivo_resumen).name}\n\n"
            f"🔬 Propiedades calculadas:\n"
            f"• Tm, ΔH, ΔG (25°C y 37°C)\n"
            f"• Contenido GC y estabilidad")
        
    except Exception as e:
        mensaje_error = f"❌ Error durante el procesamiento: {str(e)}"
        logger.error(mensaje_error)
        logger.error(traceback.format_exc())
        print(mensaje_error)
        InterfazUsuario.mostrar_mensaje("Error", mensaje_error, "error")


if __name__ == "__main__":
    main()
