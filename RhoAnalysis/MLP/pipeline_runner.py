#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Pipeline de ejecución parametrizado y modular para el análisis Rho.

Este script permite:
1. Ejecutar múltiples instancias de rhoHistFromTree.py con diferentes parámetros
2. Extraer automáticamente las rutas de archivos de salida
3. Actualizar archivos de configuración YAML (como OptimalVariable) con las nuevas rutas
4. Ejecutar CompareAlgs.py al finalizar

Uso típico:
    python RhoAnalysis/pipeline_runner.py -c config/pipeline/my_pipeline.yaml

Ejemplo de archivo de configuración del pipeline (pipeline_config.yaml):
    
    pipeline:
      # Configuración base compartida por todas las ejecuciones
      base_config:
        hist_config: '/nfs/cms/arqolmo/ExamplesFCCFullSim/config/histograms/rho_analysis_config.yml'
        reco_config: '/nfs/cms/arqolmo/ExamplesFCCFullSim/config/default/taurecolong_optimal.yaml'
        verbose: true
        
      # Lista de ejecuciones de rhoHistFromTree.py
      runs:
        - name: "Zqq_sample"
          decay: 2
          tree_file: '/nfs/cms/arqolmo/ExamplesFCCFullSim/Results/RhoAnalysis/Zqq_sample.../tree.root'
          ang: [0.5]
          meson_cut: [0.0, 45.57]
          lepton_cut: [0.0, 41.16]
          # Mapeo para archivos de config de plots
          output_mapping:
            - config_file: 'config/plots/Optimal_Variable/Rho Decay/OptimalVariableBK_d2 optimal.yaml'
              dataset_label: "Zqq BG"
              
        - name: "Bhabha_sample"
          decay: 2
          tree_file: '/path/to/bhabha_tree.root'
          ang: [3.06, 5.0]
          output_mapping:
            - config_file: 'config/plots/Optimal_Variable/Rho Decay/OptimalVariableBK_d2 optimal.yaml'
              dataset_label: "Bhabha BG"
      
      # Configuración para la ejecución final de CompareAlgs.py
      compare_algs:
        enabled: true
        config_file: 'config/plots/Optimal_Variable/Rho Decay/OptimalVariableBK_d2 optimal.yaml'
        output_dir: 'OptimalBKd2'
"""

import argparse
import os
import sys
import yaml
import subprocess
import shlex
import re
from pathlib import Path
from typing import Dict, List, Any, Optional, Tuple
from dataclasses import dataclass, field
import logging
from datetime import datetime
import numpy as np
from concurrent.futures import ProcessPoolExecutor, as_completed
import multiprocessing

# Configuración de logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


@dataclass
class RunConfig:
    """Configuración para una ejecución individual de rhoHistFromTree.py"""
    name: str
    decay: int
    tree_file: str
    hist_config: str
    reco_config: str
    ang: List[float] = field(default_factory=lambda: [0.0])
    meson_cut: List[float] = field(default_factory=lambda: [0.0, np.inf])
    lepton_cut: List[float] = field(default_factory=lambda: [0.0, np.inf])
    zmass_cut: List[float] = field(default_factory=lambda: [0.0, np.inf])
    verbose: bool = False
    output_mapping: List[Dict[str, str]] = field(default_factory=list)


def compute_output_path(
    tree_file: str,
    decay: int,
    reco_config: str,
    ang: List[float],
    meson_cut: List[float],
    lepton_cut: List[float],
    zmass_cut: List[float],
) -> str:
    """
    Calcula la ruta del archivo de salida siguiendo la misma lógica que rhoHistFromTree.py.
    
    La lógica es:
    - outputpath = directorio del tree_file
    - fileOutName = "Histos_" + sufijos de cortes + fileOutName del config
    
    Args:
        tree_file: Ruta del archivo ROOT de entrada con los árboles
        decay: Tipo de decay (0, 1, 2, etc.)
        reco_config: Archivo de configuración de reconstrucción
        ang: Rango de separación angular [min] o [min, max]
        meson_cut: Rango de corte en momento del mesón
        lepton_cut: Rango de corte en momento del leptón
        zmass_cut: Rango de corte en masa Z
        
    Returns:
        Ruta completa del archivo ROOT de salida esperado
    """
    # Cargar configuración para obtener el fileOutName base
    with open(reco_config, 'r') as f:
        config = yaml.safe_load(f)
    
    # Extraer parámetros del config
    general = config.get('general', {})
    cuts = config.get('cuts', {})
    
    outfile_base = general.get('outfile', 'output')
    
    # Construir sufijo del nombre de archivo basado en parámetros
    # Similar a myutils.setup_analysis_config
    dRMax = cuts.get('dRMax', 0.4)
    TauPhotonPCut = cuts.get('TauPhotonPCut', 0.35)
    TauPionPCut = cuts.get('TauPionPCut', 0)
    NeutronCut = cuts.get('NeutronCut', 3)
    generalPCut = cuts.get('generalPCut', 0.0)
    
    fileOutName_base = f"{outfile_base}decay{decay}_{dRMax}_tph{TauPhotonPCut}_tpi{TauPionPCut}_n{NeutronCut}_g{generalPCut}.root"
    
    # Construir el string de cortes adicionales
    out_histos_string = "Histos_"
    
    # Normalizar ang
    if len(ang) == 1:
        ang = [ang[0], np.inf]
    
    if ang[0] > 0:
        ang_max_str = str(ang[1]) if ang[1] != np.inf else "inf"
        out_histos_string += f"dRgt{ang[0]}_{ang_max_str}_"
    
    if meson_cut[0] > 0 or meson_cut[1] < 100:
        meson_max_str = str(meson_cut[1]) if meson_cut[1] != np.inf else "inf"
        out_histos_string += f"MesonPgt{meson_cut[0]}_lt{meson_max_str}_"
    
    if lepton_cut[0] > 0 or lepton_cut[1] < 100:
        lepton_max_str = str(lepton_cut[1]) if lepton_cut[1] != np.inf else "inf"
        out_histos_string += f"LeptonPgt{lepton_cut[0]}_lt{lepton_max_str}_"
    
    if zmass_cut[0] > 0 or zmass_cut[1] < 200:
        zmass_max_str = str(zmass_cut[1]) if zmass_cut[1] != np.inf else "inf"
        out_histos_string += f"Zmassgt{zmass_cut[0]}_lt{zmass_max_str}_"
    
    # Ruta de salida es el directorio del tree_file
    outputpath = os.path.dirname(tree_file)
    
    # Nombre final del archivo
    fileOutName = os.path.join(outputpath, out_histos_string + fileOutName_base)
    
    return fileOutName


def build_rho_hist_command(run_config: RunConfig) -> List[str]:
    """
    Construye el comando para ejecutar rhoHistFromTree.py.
    
    Args:
        run_config: Configuración de la ejecución
        
    Returns:
        Lista de argumentos del comando
    """
    script_path = os.path.join(
        os.path.dirname(os.path.abspath(__file__)),
        'rhoHistFromTree.py'
    )
    
    cmd = [
        'python', script_path,
        '-d', str(run_config.decay),
        '--hist-config', run_config.hist_config,
        '--tree-file', run_config.tree_file,
        '-c', run_config.reco_config,
    ]
    
    # Añadir cortes opcionales
    if run_config.ang:
        cmd.extend(['--ang'] + [str(x) for x in run_config.ang])
    
    if run_config.meson_cut and run_config.meson_cut != [0.0, np.inf]:
        cmd.extend(['--meson-cut'] + [str(x) for x in run_config.meson_cut])
    
    if run_config.lepton_cut and run_config.lepton_cut != [0.0, np.inf]:
        cmd.extend(['--lepton-cut'] + [str(x) for x in run_config.lepton_cut])
    
    if run_config.zmass_cut and run_config.zmass_cut != [0.0, np.inf]:
        cmd.extend(['--zmass-cut'] + [str(x) for x in run_config.zmass_cut])
    
    if run_config.verbose:
        cmd.append('-v')
    
    return cmd


def run_rho_hist_from_tree(run_config: RunConfig, dry_run: bool = False, quiet: bool = False) -> Tuple[bool, str, str]:
    """
    Ejecuta rhoHistFromTree.py con la configuración dada.
    
    Args:
        run_config: Configuración de la ejecución
        dry_run: Si es True, solo muestra el comando sin ejecutar
        quiet: Si es True, no muestra logs detallados (para ejecución en paralelo)
        
    Returns:
        Tuple (éxito, ruta_salida, nombre_run)
    """
    cmd = build_rho_hist_command(run_config)
    cmd_str = ' '.join(shlex.quote(c) for c in cmd)
    
    if not quiet:
        logger.info(f"[{run_config.name}] Ejecutando comando:")
        logger.info(f"  {cmd_str}")
    
    # Calcular ruta de salida esperada
    output_path = compute_output_path(
        tree_file=run_config.tree_file,
        decay=run_config.decay,
        reco_config=run_config.reco_config,
        ang=run_config.ang,
        meson_cut=run_config.meson_cut,
        lepton_cut=run_config.lepton_cut,
        zmass_cut=run_config.zmass_cut,
    )
    
    if not quiet:
        logger.info(f"[{run_config.name}] Archivo de salida esperado:")
        logger.info(f"  {output_path}")
    
    if dry_run:
        if not quiet:
            logger.info(f"[{run_config.name}] (dry-run) Comando no ejecutado")
        return True, output_path, run_config.name
    
    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            cwd=os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        )
        
        if result.returncode != 0:
            if not quiet:
                logger.error(f"[{run_config.name}] Error en ejecución:")
                logger.error(result.stderr)
            return False, output_path, run_config.name
        
        if not quiet:
            logger.info(f"[{run_config.name}] Ejecución completada exitosamente")
        
            # Verificar que el archivo de salida existe
            if os.path.exists(output_path):
                logger.info(f"[{run_config.name}] Archivo de salida verificado: {output_path}")
            else:
                logger.warning(f"[{run_config.name}] Archivo de salida no encontrado: {output_path}")
        
        return True, output_path, run_config.name
        
    except Exception as e:
        if not quiet:
            logger.error(f"[{run_config.name}] Excepción durante ejecución: {e}")
        return False, output_path, run_config.name


def _run_single_task(args: Tuple[RunConfig, bool]) -> Tuple[bool, str, str]:
    """
    Wrapper para ejecutar una tarea individual en un proceso separado.
    
    Args:
        args: Tupla (run_config, dry_run)
        
    Returns:
        Tuple (éxito, ruta_salida, nombre_run)
    """
    run_config, dry_run = args
    return run_rho_hist_from_tree(run_config, dry_run=dry_run, quiet=True)


def run_rho_hist_parallel(
    run_configs: List[RunConfig],
    dry_run: bool = False,
    max_workers: int = None,
    verbose: bool = False,
) -> Dict[str, Dict[str, Any]]:
    """
    Ejecuta múltiples instancias de rhoHistFromTree.py en paralelo.
    
    Args:
        run_configs: Lista de configuraciones de ejecución
        dry_run: Si es True, solo muestra los comandos sin ejecutar
        max_workers: Número máximo de workers (por defecto: número de CPUs)
        verbose: Si es True, muestra información adicional de progreso
        
    Returns:
        Diccionario con los resultados de cada ejecución
    """
    if max_workers is None:
        max_workers = min(len(run_configs), multiprocessing.cpu_count())
    
    logger.info(f"Ejecutando {len(run_configs)} tareas en paralelo con {max_workers} workers...")
    
    # Mostrar resumen de tareas
    if verbose:
        logger.info("Tareas a ejecutar:")
        for rc in run_configs:
            logger.info(f"  - {rc.name}: {os.path.basename(rc.tree_file)}")
    
    results = {}
    tasks = [(rc, dry_run) for rc in run_configs]
    
    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        # Enviar todas las tareas
        future_to_config = {
            executor.submit(_run_single_task, task): task[0] 
            for task in tasks
        }
        
        # Recoger resultados a medida que completan
        completed = 0
        total = len(future_to_config)
        
        for future in as_completed(future_to_config):
            run_config = future_to_config[future]
            completed += 1
            
            try:
                success, output_path, name = future.result()
                status = "✓" if success else "✗"
                logger.info(f"  [{completed}/{total}] {status} {name} completado")
                
                if verbose and not success:
                    logger.warning(f"    Falló: {output_path}")
                
                results[name] = {
                    'success': success,
                    'output_path': output_path,
                    'config': run_config,
                }
            except Exception as e:
                logger.error(f"  [{completed}/{total}] ✗ {run_config.name} falló con excepción: {e}")
                results[run_config.name] = {
                    'success': False,
                    'output_path': None,
                    'config': run_config,
                }
    
    return results


def update_yaml_config_path(
    config_file: str,
    dataset_label: str,
    new_path: str,
    base_path: str = None,
) -> bool:
    """
    Actualiza la ruta de un dataset específico en un archivo YAML de configuración.
    
    Args:
        config_file: Ruta al archivo YAML de configuración
        dataset_label: Etiqueta del dataset a actualizar
        new_path: Nueva ruta del archivo ROOT
        base_path: Ruta base para convertir a ruta relativa (opcional)
        
    Returns:
        True si la actualización fue exitosa
    """
    try:
        with open(config_file, 'r') as f:
            config = yaml.safe_load(f)
        
        if 'datasets' not in config:
            logger.error(f"No se encontró 'datasets' en {config_file}")
            return False
        
        # Convertir a ruta relativa si se proporciona base_path
        if base_path:
            try:
                new_path = os.path.relpath(new_path, base_path)
            except ValueError:
                pass  # Si no se puede hacer relativa, usar la absoluta
        
        # Buscar y actualizar el dataset con la etiqueta dada
        updated = False
        for dataset in config['datasets']:
            if dataset.get('label') == dataset_label:
                old_path = dataset.get('path', 'N/A')
                dataset['path'] = new_path
                logger.info(f"Actualizado '{dataset_label}':")
                logger.info(f"  Anterior: {old_path}")
                logger.info(f"  Nuevo:    {new_path}")
                updated = True
                # No break, por si hay múltiples datasets con la misma etiqueta
        
        if not updated:
            logger.warning(f"No se encontró dataset con label '{dataset_label}' en {config_file}")
            return False
        
        # Escribir el archivo actualizado
        with open(config_file, 'w') as f:
            yaml.dump(config, f, default_flow_style=False, allow_unicode=True, sort_keys=False)
        
        logger.info(f"Archivo de configuración actualizado: {config_file}")
        return True
        
    except Exception as e:
        logger.error(f"Error actualizando {config_file}: {e}")
        return False


def run_compare_algs(
    config_file: str,
    output_dir: str,
    dry_run: bool = False,
) -> bool:
    """
    Ejecuta TauAnalysis/CompareAlgs.py.
    
    Args:
        config_file: Archivo de configuración YAML
        output_dir: Directorio de salida para las imágenes
        dry_run: Si es True, solo muestra el comando sin ejecutar
        
    Returns:
        True si la ejecución fue exitosa
    """
    script_path = os.path.join(
        os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
        'TauAnalysis',
        'CompareAlgs.py'
    )
    
    cmd = [
        'python', script_path,
        '-c', config_file,
        '-o', output_dir,
    ]
    
    cmd_str = ' '.join(shlex.quote(c) for c in cmd)
    
    logger.info("Ejecutando CompareAlgs.py:")
    logger.info(f"  {cmd_str}")
    
    if dry_run:
        logger.info("(dry-run) Comando no ejecutado")
        return True
    
    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            cwd=os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        )
        
        if result.returncode != 0:
            logger.error("Error en CompareAlgs.py:")
            logger.error(result.stderr)
            return False
        
        logger.info("CompareAlgs.py ejecutado exitosamente")
        logger.info(result.stdout)
        return True
        
    except Exception as e:
        logger.error(f"Excepción durante CompareAlgs.py: {e}")
        return False


def load_pipeline_config(config_file: str) -> Dict[str, Any]:
    """
    Carga y valida el archivo de configuración del pipeline.
    
    Args:
        config_file: Ruta al archivo YAML del pipeline
        
    Returns:
        Diccionario con la configuración del pipeline
    """
    with open(config_file, 'r') as f:
        config = yaml.safe_load(f)
    
    if 'pipeline' not in config:
        raise ValueError("El archivo de configuración debe contener una clave 'pipeline'")
    
    return config['pipeline']


def parse_run_config(
    run_dict: Dict[str, Any],
    base_config: Dict[str, Any],
) -> RunConfig:
    """
    Parsea un diccionario de configuración de ejecución y lo combina con la configuración base.
    
    Args:
        run_dict: Diccionario con la configuración específica de la ejecución
        base_config: Configuración base compartida
        
    Returns:
        RunConfig con la configuración combinada
    """
    # Combinar configuración base con la específica
    combined = {**base_config, **run_dict}
    
    # Manejar valores por defecto para listas
    ang = combined.get('ang', [0.0])
    if isinstance(ang, (int, float)):
        ang = [ang]
    
    meson_cut = combined.get('meson_cut', [0.0, float('inf')])
    if isinstance(meson_cut, (int, float)):
        meson_cut = [meson_cut, float('inf')]
    
    lepton_cut = combined.get('lepton_cut', [0.0, float('inf')])
    if isinstance(lepton_cut, (int, float)):
        lepton_cut = [lepton_cut, float('inf')]
    
    zmass_cut = combined.get('zmass_cut', [0.0, float('inf')])
    if isinstance(zmass_cut, (int, float)):
        zmass_cut = [zmass_cut, float('inf')]
    
    return RunConfig(
        name=combined.get('name', 'unnamed'),
        decay=combined.get('decay', 2),
        tree_file=combined['tree_file'],
        hist_config=combined['hist_config'],
        reco_config=combined['reco_config'],
        ang=ang,
        meson_cut=meson_cut,
        lepton_cut=lepton_cut,
        zmass_cut=zmass_cut,
        verbose=combined.get('verbose', False),
        output_mapping=combined.get('output_mapping', []),
    )


def main():
    parser = argparse.ArgumentParser(
        description="Pipeline de ejecución parametrizado para análisis Rho",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__
    )
    
    parser.add_argument(
        '-c', '--config',
        type=str,
        required=True,
        help='Archivo YAML con la configuración del pipeline'
    )
    
    parser.add_argument(
        '--dry-run',
        action='store_true',
        help='Muestra los comandos sin ejecutarlos'
    )
    
    parser.add_argument(
        '--skip-hist',
        action='store_true',
        help='Salta la generación de histogramas (solo actualiza configs y ejecuta compare)'
    )
    
    parser.add_argument(
        '--skip-compare',
        action='store_true',
        help='Salta la ejecución de CompareAlgs.py'
    )
    
    parser.add_argument(
        '--only-run',
        type=str,
        nargs='+',
        help='Ejecuta solo las runs con los nombres especificados'
    )
    
    parser.add_argument(
        '-p', '--parallel',
        action='store_true',
        help='Ejecuta las tareas de histogramas en paralelo'
    )
    
    parser.add_argument(
        '-j', '--jobs',
        type=int,
        default=None,
        help='Número de workers para ejecución paralela (por defecto: número de CPUs)'
    )
    
    parser.add_argument(
        '-v', '--verbose',
        action='store_true',
        help='Modo verbose (muestra información adicional de progreso en paralelo)'
    )
    
    args = parser.parse_args()
    
    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)
    
    # Cargar configuración del pipeline
    logger.info(f"Cargando configuración del pipeline: {args.config}")
    try:
        pipeline_config = load_pipeline_config(args.config)
    except Exception as e:
        logger.error(f"Error cargando configuración: {e}")
        sys.exit(1)
    
    base_config = pipeline_config.get('base_config', {})
    runs = pipeline_config.get('runs', [])
    compare_config = pipeline_config.get('compare_algs', {})
    
    logger.info(f"Pipeline con {len(runs)} ejecuciones configuradas")
    
    # Preparar configuraciones de ejecución
    run_configs_to_execute = []
    
    for run_dict in runs:
        run_config = parse_run_config(run_dict, base_config)
        
        # Filtrar por nombre si se especificó --only-run
        if args.only_run and run_config.name not in args.only_run:
            logger.info(f"Saltando '{run_config.name}' (no está en --only-run)")
            continue
        
        run_configs_to_execute.append(run_config)
    
    # Ejecutar rhoHistFromTree.py (secuencial o paralelo)
    results = {}
    
    if not args.skip_hist:
        if args.parallel and len(run_configs_to_execute) > 1:
            # Ejecución en paralelo
            logger.info("=" * 60)
            logger.info("EJECUCIÓN EN PARALELO")
            logger.info("=" * 60)
            
            results = run_rho_hist_parallel(
                run_configs=run_configs_to_execute,
                dry_run=args.dry_run,
                max_workers=args.jobs,
                verbose=args.verbose,
            )
        else:
            # Ejecución secuencial (muestra toda la salida)
            for run_config in run_configs_to_execute:
                logger.info("=" * 60)
                logger.info(f"Procesando: {run_config.name}")
                logger.info("=" * 60)
                
                success, output_path, _ = run_rho_hist_from_tree(
                    run_config,
                    dry_run=args.dry_run,
                    quiet=False  # Muestra toda la salida en modo secuencial
                )
                
                results[run_config.name] = {
                    'success': success,
                    'output_path': output_path,
                    'config': run_config,
                }
    else:
        # Solo calcular rutas de salida sin ejecutar
        for run_config in run_configs_to_execute:
            output_path = compute_output_path(
                tree_file=run_config.tree_file,
                decay=run_config.decay,
                reco_config=run_config.reco_config,
                ang=run_config.ang,
                meson_cut=run_config.meson_cut,
                lepton_cut=run_config.lepton_cut,
                zmass_cut=run_config.zmass_cut,
            )
            logger.info(f"[{run_config.name}] Saltando generación de histogramas")
            logger.info(f"[{run_config.name}] Ruta de salida calculada: {output_path}")
            
            results[run_config.name] = {
                'success': True,
                'output_path': output_path,
                'config': run_config,
            }
    
    # Actualizar archivos de configuración de plots
    for name, result in results.items():
        if result['success'] and result['config'].output_mapping:
            base_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
            
            for mapping in result['config'].output_mapping:
                config_file = mapping.get('config_file')
                dataset_label = mapping.get('dataset_label')
                
                if config_file and dataset_label:
                    # Convertir ruta relativa a absoluta
                    if not os.path.isabs(config_file):
                        config_file = os.path.join(base_path, config_file)
                    
                    update_yaml_config_path(
                        config_file=config_file,
                        dataset_label=dataset_label,
                        new_path=result['output_path'],
                        base_path=base_path,
                    )
    
    # Resumen de resultados
    logger.info("=" * 60)
    logger.info("RESUMEN DE EJECUCIONES")
    logger.info("=" * 60)
    
    for name, result in results.items():
        status = "✓" if result['success'] else "✗"
        logger.info(f"  {status} {name}: {result['output_path']}")
    
    # Ejecutar CompareAlgs.py si está habilitado
    if compare_config.get('enabled', False) and not args.skip_compare:
        logger.info("=" * 60)
        logger.info("Ejecutando CompareAlgs.py")
        logger.info("=" * 60)
        
        config_file = compare_config.get('config_file')
        output_dir = compare_config.get('output_dir', 'CompareOutput')
        
        # Convertir ruta relativa a absoluta
        base_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        if config_file and not os.path.isabs(config_file):
            config_file = os.path.join(base_path, config_file)
        
        run_compare_algs(
            config_file=config_file,
            output_dir=output_dir,
            dry_run=args.dry_run,
        )
    
    logger.info("Pipeline completado")


if __name__ == '__main__':
    main()
