#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Visualización de Resultados - Método de Momentos
Parche Conductor en Interfaz Dieléctrica

Este script genera visualizaciones 2D y 3D de:
- Distribución de carga superficial
- Campo de potencial eléctrico
- Campo eléctrico (magnitud y vectores)
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
import seaborn as sns

# Configuración de estilo
plt.style.use('seaborn-v0_8-darkgrid')
sns.set_palette("husl")
plt.rcParams['figure.figsize'] = (12, 8)
plt.rcParams['font.size'] = 11

# ==============================================================================
# CARGAR DATOS
# ==============================================================================

print("CARGA DE DATOS")
charge_df = pd.read_csv('charge_distribution.csv')
field_df = pd.read_csv('field_data.csv')

print(f"Distribución de carga: {len(charge_df)} puntos")
print(f"Datos de campo: {len(field_df)} puntos")

# ==============================================================================
# 1. DISTRIBUCIÓN DE CARGA SUPERFICIAL
# ==============================================================================

def plot_charge_distribution():
    """Gráfica 2D y 3D de la distribución de carga"""
    
    fig = plt.figure(figsize=(20, 8))
    
    # --- Subplot 1: Vista 2D (mapa de calor) ---
    ax1 = fig.add_subplot(131)
    
    # Crear malla para contour plot
    x_unique = np.sort(charge_df['x'].unique())
    y_unique = np.sort(charge_df['y'].unique())
    X, Y = np.meshgrid(x_unique, y_unique)
    
    # Reorganizar sigma en matriz
    Z = charge_df.pivot(index='y', columns='x', values='sigma').values
    
    # Mapa de calor
    im = ax1.contourf(X * 100, Y * 100, Z * 1e9, levels=20, cmap='YlOrRd')
    ax1.contour(X * 100, Y * 100, Z * 1e9, levels=10, colors='black', 
                linewidths=0.5, alpha=0.3)
    
    cbar = plt.colorbar(im, ax=ax1)
    cbar.set_label('Densidad de Carga σ (nC/m²)', fontsize=12)
    
    ax1.set_xlabel('x (cm)', fontsize=12)
    ax1.set_ylabel('y (cm)', fontsize=12)
    ax1.set_title('Distribución de Carga Superficial\n(Vista 2D)', fontsize=14, fontweight='bold')
    ax1.set_aspect('equal')
    ax1.grid(True, alpha=0.3)
    
    # --- Subplot 2: Vista 3D ---
    ax2 = fig.add_subplot(132, projection='3d')
    
    surf = ax2.plot_surface(X * 100, Y * 100, Z * 1e9, cmap='YlOrRd', 
                            edgecolor='none', alpha=0.9)
    
    ax2.set_xlabel('x (cm)', fontsize=11)
    ax2.set_ylabel('y (cm)', fontsize=11)
    ax2.set_zlabel('σ (nC/m²)', fontsize=11)
    ax2.set_title('Distribución de Carga\n(Vista 3D)', fontsize=14, fontweight='bold')
    ax2.view_init(elev=25, azim=45)
    
    # --- Subplot 3: Cortes a lo largo de líneas centrales ---
    ax3 = fig.add_subplot(133)
    
    # Corte en y=0 (línea central horizontal)
    center_y_idx = len(y_unique) // 2
    charge_x = charge_df[np.isclose(charge_df['y'], y_unique[center_y_idx], atol=1e-6)]
    charge_x = charge_x.sort_values('x')
    
    # Corte en x=0 (línea central vertical)
    center_x_idx = len(x_unique) // 2
    charge_y = charge_df[np.isclose(charge_df['x'], x_unique[center_x_idx], atol=1e-6)]
    charge_y = charge_y.sort_values('y')
    
    ax3.plot(charge_x['x'] * 100, charge_x['sigma'] * 1e9, 'o-', 
             linewidth=2, markersize=6, label='Corte en y=0')
    ax3.plot(charge_y['y'] * 100, charge_y['sigma'] * 1e9, 's-', 
             linewidth=2, markersize=6, label='Corte en x=0')
    
    ax3.set_xlabel('Posición (cm)', fontsize=12)
    ax3.set_ylabel('Densidad de Carga σ (nC/m²)', fontsize=12)
    ax3.set_title('Cortes por Líneas Centrales', fontsize=14, fontweight='bold')
    ax3.legend(fontsize=11)
    ax3.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('distribucion_carga.png', dpi=300, bbox_inches='tight')
    print("Guardado: distribucion_carga.png")
    plt.show()

# ==============================================================================
# 2. ESTADÍSTICAS DE CARGA
# ==============================================================================

def plot_charge_statistics():
    """Histograma y estadísticas de la distribución de carga"""
    
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    
    # --- Subplot 1: Histograma ---
    ax1 = axes[0]
    
    n, bins, patches = ax1.hist(charge_df['sigma'] * 1e9, bins=30, 
                                 color='steelblue', edgecolor='black', alpha=0.7)
    
    # Colorear barras según valor
    cm_map = cm.get_cmap('YlOrRd')
    bin_centers = 0.5 * (bins[:-1] + bins[1:])
    col = bin_centers - min(bin_centers)
    col /= max(col)
    
    for c, p in zip(col, patches):
        plt.setp(p, 'facecolor', cm_map(c))
    
    ax1.axvline(charge_df['sigma'].mean() * 1e9, color='red', 
                linestyle='--', linewidth=2, label='Media')
    ax1.axvline(charge_df['sigma'].median() * 1e9, color='green', 
                linestyle='--', linewidth=2, label='Mediana')
    
    ax1.set_xlabel('Densidad de Carga σ (nC/m²)', fontsize=12)
    ax1.set_ylabel('Frecuencia', fontsize=12)
    ax1.set_title('Histograma de Densidad de Carga', fontsize=14, fontweight='bold')
    ax1.legend(fontsize=11)
    ax1.grid(True, alpha=0.3)
    
    # --- Subplot 2: Box plot ---
    ax2 = axes[1]
    
    bp = ax2.boxplot([charge_df['sigma'] * 1e9], 
                      vert=True, patch_artist=True, widths=0.5)
    
    bp['boxes'][0].set_facecolor('lightblue')
    bp['boxes'][0].set_edgecolor('black')
    bp['medians'][0].set_color('red')
    bp['medians'][0].set_linewidth(2)
    
    ax2.set_ylabel('Densidad de Carga σ (nC/m²)', fontsize=12)
    ax2.set_title('Diagrama de Caja', fontsize=14, fontweight='bold')
    ax2.grid(True, alpha=0.3, axis='y')
    ax2.set_xticklabels(['Distribución'])
    
    # Añadir estadísticas
    stats_text = f"""Estadísticas:
    Media: {charge_df['sigma'].mean()*1e9:.3f} nC/m²
    Mediana: {charge_df['sigma'].median()*1e9:.3f} nC/m²
    Desv. Est.: {charge_df['sigma'].std()*1e9:.3f} nC/m²
    Mín: {charge_df['sigma'].min()*1e9:.3f} nC/m²
    Máx: {charge_df['sigma'].max()*1e9:.3f} nC/m²
    """
    
    ax2.text(1.5, charge_df['sigma'].mean() * 1e9, stats_text,
             bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5),
             fontsize=10, verticalalignment='center')
    
    plt.tight_layout()
    plt.savefig('estadisticas_carga.png', dpi=300, bbox_inches='tight')
    print("Guardado: estadisticas_carga.png")
    plt.show()

# ==============================================================================
# 3. CAMPO DE POTENCIAL ELÉCTRICO
# ==============================================================================

def plot_potential_field():
    """Visualización del potencial eléctrico en varios planos"""
    
    # Seleccionar planos z específicos
    z_planes = [0.01, 0.03, 0.05]  # 1cm, 3cm, 5cm sobre la interfaz
    
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))
    
    for idx, z_val in enumerate(z_planes):
        ax = axes[idx]
        
        # Filtrar datos para este plano z
        plane_data = field_df[np.isclose(field_df['z'], z_val, atol=0.001)]
        
        if len(plane_data) > 0:
            # Crear malla
            x_unique = np.sort(plane_data['x'].unique())
            y_unique = np.sort(plane_data['y'].unique())
            
            if len(x_unique) > 1 and len(y_unique) > 1:
                X, Y = np.meshgrid(x_unique, y_unique)
                V_matrix = plane_data.pivot_table(
                    index='y', columns='x', values='V', aggfunc='mean'
                ).values
                
                # Mapa de calor
                im = ax.contourf(X * 100, Y * 100, V_matrix, levels=15, cmap='coolwarm')
                ax.contour(X * 100, Y * 100, V_matrix, levels=10, 
                          colors='black', linewidths=0.5, alpha=0.4)
                
                cbar = plt.colorbar(im, ax=ax)
                cbar.set_label('Potencial V (V)', fontsize=10)
                
                ax.set_xlabel('x (cm)', fontsize=11)
                ax.set_ylabel('y (cm)', fontsize=11)
                ax.set_title(f'Potencial en z = {z_val*100:.1f} cm', 
                            fontsize=12, fontweight='bold')
                ax.set_aspect('equal')
                ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('campo_potencial.png', dpi=300, bbox_inches='tight')
    print("Guardado: campo_potencial.png")
    plt.show()

# ==============================================================================
# 4. POTENCIAL EN PLANO XZ (CORTE TRANSVERSAL)
# ==============================================================================

def plot_potential_xz_plane():
    """Potencial en plano xz (corte y=0)"""
    
    # Filtrar datos cerca de y=0
    xz_data = field_df[np.abs(field_df['y']) < 0.005]
    
    if len(xz_data) > 0:
        fig = plt.figure(figsize=(14, 6))
        
        # --- Subplot 1: Contour plot ---
        ax1 = fig.add_subplot(121)
        
        x_unique = np.sort(xz_data['x'].unique())
        z_unique = np.sort(xz_data['z'].unique())
        X, Z = np.meshgrid(x_unique, z_unique)
        
        V_matrix = xz_data.pivot_table(
            index='z', columns='x', values='V', aggfunc='mean'
        ).values
        
        im = ax1.contourf(X * 100, Z * 100, V_matrix, levels=20, cmap='viridis')
        ax1.contour(X * 100, Z * 100, V_matrix, levels=10, 
                   colors='white', linewidths=0.5, alpha=0.5)
        
        # Marcar la interfaz
        ax1.axhline(y=0, color='red', linestyle='--', linewidth=2, 
                   label='Interfaz Dieléctrica')
        
        cbar = plt.colorbar(im, ax=ax1)
        cbar.set_label('Potencial V (V)', fontsize=12)
        
        ax1.set_xlabel('x (cm)', fontsize=12)
        ax1.set_ylabel('z (cm)', fontsize=12)
        ax1.set_title('Potencial en Plano xz (y=0)', fontsize=14, fontweight='bold')
        ax1.legend(fontsize=10)
        ax1.grid(True, alpha=0.3)
        
        # --- Subplot 2: Corte vertical en x=0 ---
        ax2 = fig.add_subplot(122)
        
        vertical_cut = xz_data[np.abs(xz_data['x']) < 0.005].sort_values('z')
        
        if len(vertical_cut) > 0:
            ax2.plot(vertical_cut['V'], vertical_cut['z'] * 100, 'o-', 
                    linewidth=2, markersize=6, color='darkblue')
            
            ax2.axhline(y=0, color='red', linestyle='--', linewidth=2, 
                       label='Interfaz (z=0)')
            ax2.axhline(y=0, color='red', linestyle='--', linewidth=2)
            
            ax2.set_xlabel('Potencial V (V)', fontsize=12)
            ax2.set_ylabel('z (cm)', fontsize=12)
            ax2.set_title('Potencial a lo Largo del Eje z\n(x=0, y=0)', 
                         fontsize=14, fontweight='bold')
            ax2.legend(fontsize=10)
            ax2.grid(True, alpha=0.3)
        
        plt.tight_layout()
        plt.savefig('potencial_xz.png', dpi=300, bbox_inches='tight')
        print("Guardado: potencial_xz.png")
        plt.show()

# ==============================================================================
# 5. CAMPO ELÉCTRICO - MAGNITUD
# ==============================================================================

def plot_electric_field_magnitude():
    """Magnitud del campo eléctrico en varios planos"""
    
    z_planes = [0.01, 0.03, 0.05]  # 1cm, 3cm, 5cm
    
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))
    
    for idx, z_val in enumerate(z_planes):
        ax = axes[idx]
        
        plane_data = field_df[np.isclose(field_df['z'], z_val, atol=0.001)]
        
        if len(plane_data) > 0:
            x_unique = np.sort(plane_data['x'].unique())
            y_unique = np.sort(plane_data['y'].unique())
            
            if len(x_unique) > 1 and len(y_unique) > 1:
                X, Y = np.meshgrid(x_unique, y_unique)
                E_matrix = plane_data.pivot_table(
                    index='y', columns='x', values='E_mag', aggfunc='mean'
                ).values
                
                # Usar escala logarítmica para mejor visualización
                im = ax.contourf(X * 100, Y * 100, np.log10(E_matrix + 1e-10), 
                                levels=15, cmap='hot')
                
                cbar = plt.colorbar(im, ax=ax)
                cbar.set_label('log₁₀(|E|) (V/m)', fontsize=10)
                
                ax.set_xlabel('x (cm)', fontsize=11)
                ax.set_ylabel('y (cm)', fontsize=11)
                ax.set_title(f'Magnitud de Campo E en z = {z_val*100:.1f} cm', 
                            fontsize=12, fontweight='bold')
                ax.set_aspect('equal')
                ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('campo_electrico_magnitud.png', dpi=300, bbox_inches='tight')
    print("Guardado: campo_electrico_magnitud.png")
    plt.show()

# ==============================================================================
# 6. VECTORES DE CAMPO ELÉCTRICO
# ==============================================================================

def plot_electric_field_vectors():
    """Vectores del campo eléctrico en plano xy"""
    
    z_val = 0.02  # 2 cm sobre la interfaz
    
    plane_data = field_df[np.isclose(field_df['z'], z_val, atol=0.001)]
    
    # Submuestrear para mejor visualización
    step = 3
    plane_data_sub = plane_data.iloc[::step]
    
    fig = plt.figure(figsize=(12, 10))
    
    # --- Campo vectorial 2D ---
    ax1 = fig.add_subplot(211)
    
    # Magnitud de fondo
    x_unique = np.sort(plane_data['x'].unique())
    y_unique = np.sort(plane_data['y'].unique())
    X, Y = np.meshgrid(x_unique, y_unique)
    E_matrix = plane_data.pivot_table(
        index='y', columns='x', values='E_mag', aggfunc='mean'
    ).values
    
    im = ax1.contourf(X * 100, Y * 100, E_matrix, levels=15, cmap='Blues', alpha=0.6)
    
    # Vectores
    Q = ax1.quiver(plane_data_sub['x'] * 100, plane_data_sub['y'] * 100,
                   plane_data_sub['Ex'], plane_data_sub['Ey'],
                   plane_data_sub['E_mag'],
                   cmap='autumn', scale=50, width=0.003, alpha=0.8)
    
    ax1.set_xlabel('x (cm)', fontsize=12)
    ax1.set_ylabel('y (cm)', fontsize=12)
    ax1.set_title(f'Campo Eléctrico Vectorial en z = {z_val*100:.1f} cm', 
                 fontsize=14, fontweight='bold')
    ax1.set_aspect('equal')
    ax1.grid(True, alpha=0.3)
    
    cbar = plt.colorbar(Q, ax=ax1)
    cbar.set_label('|E| (V/m)', fontsize=11)
    
    # --- Campo vectorial en plano xz (y=0) ---
    ax2 = fig.add_subplot(212)
    
    xz_data = field_df[np.abs(field_df['y']) < 0.005]
    xz_sub = xz_data.iloc[::step*2]
    
    # Magnitud de fondo
    x_unique_xz = np.sort(xz_data['x'].unique())
    z_unique_xz = np.sort(xz_data['z'].unique())
    X_xz, Z_xz = np.meshgrid(x_unique_xz, z_unique_xz)
    E_matrix_xz = xz_data.pivot_table(
        index='z', columns='x', values='E_mag', aggfunc='mean'
    ).values
    
    im2 = ax2.contourf(X_xz * 100, Z_xz * 100, E_matrix_xz, 
                       levels=15, cmap='Blues', alpha=0.6)
    
    # Vectores
    Q2 = ax2.quiver(xz_sub['x'] * 100, xz_sub['z'] * 100,
                    xz_sub['Ex'], xz_sub['Ez'],
                    xz_sub['E_mag'],
                    cmap='autumn', scale=50, width=0.003, alpha=0.8)
    
    ax2.axhline(y=0, color='red', linestyle='--', linewidth=2, 
               label='Interfaz Dieléctrica')
    
    ax2.set_xlabel('x (cm)', fontsize=12)
    ax2.set_ylabel('z (cm)', fontsize=12)
    ax2.set_title('Campo Eléctrico en Plano xz (y=0)', 
                 fontsize=14, fontweight='bold')
    ax2.legend(fontsize=10)
    ax2.grid(True, alpha=0.3)
    
    cbar2 = plt.colorbar(Q2, ax=ax2)
    cbar2.set_label('|E| (V/m)', fontsize=11)
    
    plt.tight_layout()
    plt.savefig('campo_electrico_vectores.png', dpi=300, bbox_inches='tight')
    print("Guardado: campo_electrico_vectores.png")
    plt.show()

# ==============================================================================
# 7. CAMPO ELÉCTRICO 3D
# ==============================================================================

def plot_electric_field_3d():
    """Visualización 3D del campo eléctrico"""
    
    # Tomar un subconjunto de puntos
    field_sub = field_df.iloc[::20]
    
    fig = plt.figure(figsize=(14, 10))
    
    # --- Subplot 1: Magnitud en 3D ---
    ax1 = fig.add_subplot(121, projection='3d')
    
    scatter = ax1.scatter(field_sub['x'] * 100, 
                         field_sub['y'] * 100, 
                         field_sub['z'] * 100,
                         c=np.log10(field_sub['E_mag'] + 1e-10),
                         cmap='plasma', s=20, alpha=0.6)
    
    ax1.set_xlabel('x (cm)', fontsize=11)
    ax1.set_ylabel('y (cm)', fontsize=11)
    ax1.set_zlabel('z (cm)', fontsize=11)
    ax1.set_title('Magnitud del Campo Eléctrico\n(Vista 3D)', 
                 fontsize=13, fontweight='bold')
    
    cbar = plt.colorbar(scatter, ax=ax1, shrink=0.5)
    cbar.set_label('log₁₀(|E|) (V/m)', fontsize=10)
    
    # --- Subplot 2: Líneas de campo ---
    ax2 = fig.add_subplot(122, projection='3d')
    
    # Tomar aún menos puntos para líneas de campo
    field_lines = field_df.iloc[::50]
    
    ax2.quiver(field_lines['x'] * 100,
               field_lines['y'] * 100,
               field_lines['z'] * 100,
               field_lines['Ex'],
               field_lines['Ey'],
               field_lines['Ez'],
               length=0.5, normalize=True, 
               color='steelblue', alpha=0.6, linewidth=0.8)
    
    ax2.set_xlabel('x (cm)', fontsize=11)
    ax2.set_ylabel('y (cm)', fontsize=11)
    ax2.set_zlabel('z (cm)', fontsize=11)
    ax2.set_title('Líneas de Campo Eléctrico\n(Vista 3D)', 
                 fontsize=13, fontweight='bold')
    
    plt.tight_layout()
    plt.savefig('campo_electrico_3d.png', dpi=300, bbox_inches='tight')
    print("Guardado: campo_electrico_3d.png")
    plt.show()

# ==============================================================================
# EJECUTAR TODAS LAS VISUALIZACIONES
# ==============================================================================

if __name__ == "__main__":
    print("\n" + "="*70)
    print("VISUALIZACIÓN DE RESULTADOS - MÉTODO DE MOMENTOS")
    print("="*70 + "\n")
    
    print("Generando visualizaciones...\n")
    
    plot_charge_distribution()
    
    plot_charge_statistics()
    
    plot_potential_field()
    
    plot_potential_xz_plane()
    
    plot_electric_field_magnitude()
    
    plot_electric_field_vectors()
    
    plot_electric_field_3d()
    
    print("="*70)
    print("\nArchivos generados:")
    print("  1. distribucion_carga.png")
    print("  2. estadisticas_carga.png")
    print("  3. campo_potencial.png")
    print("  4. potencial_xz.png")
    print("  5. campo_electrico_magnitud.png")
    print("  6. campo_electrico_vectores.png")
    print("  7. campo_electrico_3d.png")
    print("\n")
