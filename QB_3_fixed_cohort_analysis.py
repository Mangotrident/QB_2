"""
QB_3: Fixed Cohort Analysis
Proper cohort statistics with real data, p-values, and effect sizes.
"""

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu
from typing import Dict, List, Optional, Tuple
from pathlib import Path
import matplotlib.pyplot as plt
import seaborn as sns


class FixedCohortAnalyzer:
    """
    Fixed cohort analyzer with proper statistics.
    """
    
    def __init__(self):
        self.results = {}
    
    def compare_groups(self, results_df: pd.DataFrame, group_column: str,
                      metric_column: str, group1_label: str, group2_label: str) -> Dict:
        """
        Compare two groups with proper statistics.
        
        Returns: p-value, Cohen's d, direction, etc.
        """
        group1 = results_df[results_df[group_column] == group1_label][metric_column].dropna()
        group2 = results_df[results_df[group_column] == group2_label][metric_column].dropna()
        
        if len(group1) == 0 or len(group2) == 0:
            return {'error': 'Insufficient data', 'n1': len(group1), 'n2': len(group2)}
        
        # Mann-Whitney U test (two-sided)
        statistic, p_value = mannwhitneyu(group1, group2, alternative='two-sided')
        
        # Cohen's d
        pooled_std = np.sqrt(
            ((len(group1) - 1) * group1.std()**2 + 
             (len(group2) - 1) * group2.std()**2) / 
            (len(group1) + len(group2) - 2)
        )
        if pooled_std > 0:
            cohens_d = (group1.mean() - group2.mean()) / pooled_std
        else:
            cohens_d = 0.0
        
        # Direction
        direction = 'Higher' if group1.mean() > group2.mean() else 'Lower'
        direction_match = True  # Will be set based on expected
        
        results = {
            'metric': metric_column,
            'group1_label': group1_label,
            'group2_label': group2_label,
            'n1': len(group1),
            'n2': len(group2),
            'mean1': float(group1.mean()),
            'mean2': float(group2.mean()),
            'std1': float(group1.std()),
            'std2': float(group2.std()),
            'median1': float(group1.median()),
            'median2': float(group2.median()),
            'statistic': float(statistic),
            'p_value': float(p_value),
            'cohens_d': float(cohens_d),
            'direction': direction,
            'significant': p_value < 0.05,
            'effect_size_adequate': abs(cohens_d) >= 0.6,
            'pass': p_value < 0.05 and abs(cohens_d) >= 0.6
        }
        
        return results
    
    def plot_comparison(self, results_df: pd.DataFrame, group_column: str,
                       metric_column: str, output_path: Path,
                       group1_label: str = 'Low', group2_label: str = 'High'):
        """
        Plot comparison with p-value and Cohen's d annotated.
        """
        fig, ax = plt.subplots(figsize=(8, 6))
        
        group1 = results_df[results_df[group_column] == group1_label][metric_column].dropna()
        group2 = results_df[results_df[group_column] == group2_label][metric_column].dropna()
        
        # Violin plot
        data_to_plot = [group1.values, group2.values]
        parts = ax.violinplot(data_to_plot, positions=[0, 1], widths=0.6,
                             showmeans=True, showmedians=True)
        
        # Color violins
        for pc in parts['bodies']:
            pc.set_facecolor('lightblue')
            pc.set_alpha(0.7)
        
        ax.set_xticks([0, 1])
        ax.set_xticklabels([group1_label, group2_label], fontsize=12)
        ax.set_ylabel(metric_column, fontsize=12, fontweight='bold')
        ax.set_title(f'{metric_column} Comparison: {group1_label} vs {group2_label}',
                    fontsize=14, fontweight='bold')
        
        # Get statistics
        stats = self.compare_groups(results_df, group_column, metric_column,
                                   group1_label, group2_label)
        
        # Annotate p-value and Cohen's d
        p_val = stats.get('p_value', 1.0)
        cohens_d = stats.get('cohens_d', 0.0)
        
        textstr = f'p = {p_val:.4f}\nCohen\'s d = {cohens_d:.3f}'
        if stats.get('pass', False):
            textstr += '\n✓ PASS'
        else:
            textstr += '\n✗ FAIL'
        
        ax.text(0.5, 0.95, textstr, transform=ax.transAxes,
               verticalalignment='top', bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5),
               fontsize=11)
        
        ax.grid(True, alpha=0.3, axis='y')
        plt.tight_layout()
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        plt.close()
    
    def analyze_complexI_split(self, results_df: pd.DataFrame,
                              complexI_scores: pd.Series,
                              metric_column: str = 'ETE_peak') -> Dict:
        """
        Analyze Low vs High Complex I quartiles.
        """
        # Add Complex I scores
        results_df = results_df.copy()
        results_df['ComplexI_score'] = results_df['sample_id'].map(complexI_scores.to_dict())
        
        # Split by quartiles
        q1 = complexI_scores.quantile(0.25)
        q4 = complexI_scores.quantile(0.75)
        
        results_df['ComplexI_group'] = 'Middle'
        results_df.loc[results_df['ComplexI_score'] <= q1, 'ComplexI_group'] = 'Low'
        results_df.loc[results_df['ComplexI_score'] >= q4, 'ComplexI_group'] = 'High'
        
        # Compare Low vs High
        comparison = self.compare_groups(results_df, 'ComplexI_group', metric_column,
                                        'Low', 'High')
        
        # Expected: High Complex I should have higher ETE_peak
        expected_direction = 'Higher'
        direction_match = comparison.get('direction') == expected_direction
        
        comparison['expected_direction'] = expected_direction
        comparison['direction_match'] = direction_match
        comparison['pass'] = (comparison.get('pass', False) and direction_match)
        
        return comparison, results_df

