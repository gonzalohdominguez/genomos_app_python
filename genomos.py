"""
GENOMOS - Análisis de genotipos por Tm (HRM-PCR)
------------------------------------------------

Este script permite clasificar genotipos a partir de datos de melting temperature (Tm)
obtenidos por HRM-PCR. Soporta uno o varios loci/mutaciones y genera:

1. Un archivo Excel con la clasificación por muestra.
2. Opcionalmente, un archivo TXT con la distribución de genotipos o estados.

Autores: Gonzalo Hernán Domínguez (CENEXA-CREG, UNLP-CONICET)
Fecha: 2025
"""

import pandas as pd
import numpy as np
import os
import argparse
from collections import Counter

def main():
    """
    Función principal del script. Se encarga de:
    1. Parsear argumentos de línea de comando.
    2. Leer el archivo Excel con datos de Tm.
    3. Validar las mutaciones y Tm de referencia ingresados por el usuario.
    4. Clasificar cada muestra según los valores de Tm más cercanos.
    5. Generar archivo Excel con los resultados.
    6. (Opcional) Generar archivo TXT con distribución y resumen.
    """

    # ================================
    # Configuración del parser de argumentos
    # ================================
    parser = argparse.ArgumentParser(
        description=(
            "GENOMOS: análisis de genotipos por Tm.\n\n"
            "IMPORTANTE:\n"
            "1. El archivo Excel debe contener una columna por mutación "
            "con el nombre 'Tm_{mutación}', por ejemplo: 'Tm_1016', 'Tm_1534'.\n"
            "2. El valor en -n (cantidad de mutaciones) debe coincidir con el número de argumentos -t.\n"
            "3. Cada -t debe tener el formato: POS:S:VAL,H:VAL,R:VAL\n"
            "   Ejemplo: -t 1016:S:73.2,H:72.66,R:72.21\n"
            "            -t 1534:S:81.71,H:81.81,R:82.36\n"
            "   Donde 'POS' debe coincidir con el sufijo de la columna 'Tm_POS' en el Excel.\n"
            "4. Opcionalmente, usar --txt <archivo.txt> para guardar distribución de genotipos y resumen por alelo."
        ),
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument("-n", "--num_mutaciones", type=int, required=True, help="Cantidad de mutaciones a analizar.")
    parser.add_argument("-f", "--file", type=str, required=True, help="Ruta al archivo Excel con las muestras.")
    parser.add_argument("-t", "--tm", action="append", required=True, help="Tm de referencia por mutación. Formato: 'POS:S:VAL,H:VAL,R:VAL'.")
    parser.add_argument("-o", "--output", type=str, default="resultados.xlsx", help="Archivo Excel de salida.")
    parser.add_argument("--txt", type=str, help="Archivo .txt para guardar distribución y resumen por alelo")
    args = parser.parse_args()

    # ================================
    # Validaciones iniciales
    # ================================
    
    # Validar número de mutaciones ingresadas (validar cantidad de -t con -n)
    if args.num_mutaciones != len(args.tm):
        raise ValueError(f"Error: -n={args.num_mutaciones} pero se definieron {len(args.tm)} mutaciones con -t.")

    # Validar archivo de entrada
    if not os.path.exists(args.file):
        raise FileNotFoundError(f"Archivo '{args.file}' no encontrado.")

    # Leer datos desde Excel
    df = pd.read_excel(args.file)

    # ================================
    # Procesamiento de Tm de referencia
    # ================================

    posiciones = {}
    for tm_entry in args.tm:
        try:
            # Ejemplo de entrada: "1016:S:73.2,H:72.66,R:72.21"
            pos, estados_str = tm_entry.split(":", 1)
            posiciones[pos] = {}

            # Parsear estados
            estados = estados_str.split(",")
            for est in estados:
                estado, valor = est.split(":")
                estado = estado.strip()
                valor = float(valor)

                # Normalizar claves de estados
                if estado.lower().startswith("s"):
                    estado_clave = "Sensible"
                elif estado.lower().startswith("h"):
                    estado_clave = "Heterocigoto"
                elif estado.lower().startswith("r"):
                    estado_clave = "Resistente"
                else:
                    raise ValueError(f"Estado '{estado}' no reconocido.")
                posiciones[pos][estado_clave] = valor

        except Exception as e:
            raise ValueError(f"Error en formato de -t '{tm_entry}': {e}")

    # ================================
    # Verificación de columnas requeridas
    # ================================

    # Verificar columnas y mantener solo las especificadas en -t
    columnas_usadas = []
    for pos in posiciones:
        col_name = f'Tm_{pos}'
        if col_name not in df.columns:
            raise ValueError(f"No se encontró la columna '{col_name}' en el archivo.")
        df[col_name] = pd.to_numeric(df[col_name], errors='coerce') # Asegurar valores numéricos
        columnas_usadas.append(col_name)

    # Reducir df a solo columnas de interés
    df = df[columnas_usadas]

    # ================================
    # Cálculo de diferencias absolutas con Tm de referencia
    # ================================

    for pos, tm_dict in posiciones.items():
        for estado, tm_val in tm_dict.items():
            df[f'{pos}_{estado}'] = (df[f'Tm_{pos}'] - tm_val).abs()

    # ================================
    # Asignar estado según Tm más cercano
    # ================================
    
    # Asignar genotipo por posición
    def genotipo_por_posicion(row, pos):
        """Dado un registro (fila) y una posición, retorna el estado (S, H o R) 
        más cercano en Tm."""

        if pd.isna(row[f'Tm_{pos}']):
            return None
        diffs = {estado: row[f'{pos}_{estado}'] for estado in posiciones[pos].keys()}
        return min(diffs, key=diffs.get) # Estado con menor diferencia

    for pos in posiciones:
        df[f'Estado_{pos}'] = df.apply(lambda row: genotipo_por_posicion(row, pos), axis=1)

    # ================================
    # Generar genotipo resultante (si hay múltiples mutaciones)
    # ================================

    def genotipo_resultante_dinamico(row):
        """
        Construye el genotipo resultante concatenando estados por posición.
        Ejemplo: S, H1, R2 → "SH1R2"
        """
        genotipo = []
        for idx, pos in enumerate(posiciones):
            estado = row[f'Estado_{pos}']
            if estado is None:
                return "No se pudo determinar"
            if estado == "Sensible":
                genotipo.append("S")
            elif estado == "Heterocigoto":
                genotipo.append(f"H{idx+1}")
            elif estado == "Resistente":
                genotipo.append(f"R{idx+1}")
        return ''.join(genotipo)

    cols_resultado = [f'Tm_{pos}' for pos in posiciones] + [f'Estado_{pos}' for pos in posiciones]

    if len(posiciones) > 1:
        df['Genotipo_Resultante'] = df.apply(genotipo_resultante_dinamico, axis=1)
        cols_resultado.append('Genotipo_Resultante')

    # Guardar archivo Excel de resultados
    resultados_df = df[cols_resultado]
    resultados_df.to_excel(args.output, index=False)
    print(f"Resultados guardados en '{args.output}'.")

    # ================================
    # Generar archivo TXT opcional
    # ================================

    # Guardar distribución y resumen por alelo
    if args.txt:
        with open(args.txt, 'w', encoding='utf-8') as f:
            if len(posiciones) > 1:
                # Caso: múltiples mutaciones
                total_muestras = len(df)
                conteo_genotipos = Counter(df['Genotipo_Resultante'])

                # Alelos posibles: H, R numerados + S
                alelos = []
                for idx in range(len(posiciones)):
                    alelos += [f"H{idx+1}", f"R{idx+1}"]
                alelos.append("S")

                # Recuento por alelo
                recuento_alelos = {a:0 for a in alelos}
                for g in df['Genotipo_Resultante']:
                    if g != "No se pudo determinar":
                        for a in alelos:
                            recuento_alelos[a] += g.count(a)

                # Escribir distribución de genotipos
                f.write("=== Distribución de Genotipos ===\n")
                f.write("Genotipo\tCantidad\tPorcentaje\n")
                for g, c in conteo_genotipos.items():
                    f.write(f"{g}\t{c}\t{(c/total_muestras)*100:.2f}%\n")

                # Escribir resumen por alelo
                f.write("\n=== Suma y porcentaje por alelo ===\n")
                f.write("Alelo\tCantidad\tPorcentaje\n")
                total_ale = sum(recuento_alelos.values())
                for a in alelos:
                    f.write(f"{a}\t{recuento_alelos[a]}\t{(recuento_alelos[a]/total_ale*100 if total_ale>0 else 0):.2f}%\n")

                print(f"Distribución y resumen por alelo guardados en '{args.txt}'.")

            else:
                # Caso: una sola mutación → contar estados simples
                pos = list(posiciones.keys())[0]
                estados_col = f'Estado_{pos}'
                total_muestras = len(df)
                conteo_estados = Counter(df[estados_col].dropna())

                # Escribir distribución de estados
                f.write(f"=== Distribución de Estados para {pos} ===\n")
                f.write("Estado\tCantidad\tPorcentaje\n")
                for estado, c in conteo_estados.items():
                    f.write(f"{estado}\t{c}\t{(c/total_muestras)*100:.2f}%\n")

                print(f"Distribución de estados guardada en '{args.txt}'.")

if __name__ == "__main__":
    main()