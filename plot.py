import json
import os
from markupsafe import escape
import html
import re
import logging
from typing import Dict, List, Tuple, Any


class Plot:
    """
    Plotting and visualization utilities for G.O.A.T

    Handles data preparation and processing for various plot types:
    - Scatter plots (2-color, multi-dataset)
    - Gene type comparisons
    - Expression heatmaps
    - Enrichment visualizations
    """

    def __init__(self):
        print('Plot loading')
        self.logger = logging.getLogger(__name__)

    def get_data(self, formdata: Dict) -> Tuple[Dict, str]:
        """
        Main entry point for plot data retrieval

        Args:
            formdata: Form data containing plot parameters

        Returns:
            Tuple of (data_dict, message)
        """
        plot_type = formdata.get('plot_type', 'scatter')

        try:
            if plot_type == 'scatter':
                return self.prepare_scatter_data(formdata)
            elif plot_type == 'scatter_colored':
                return self.prepare_colored_scatter_data(formdata)
            elif plot_type == 'gene_type_comparison':
                return self.prepare_gene_type_comparison(formdata)
            elif plot_type == 'heatmap':
                return self.prepare_heatmap_data(formdata)
            else:
                return {}, f'Unknown plot type: {plot_type}'
        except Exception as e:
            self.logger.error(f'Error in get_data: {e}')
            return {}, f'Error preparing plot data: {str(e)}'

    def prepare_scatter_data(self, formdata: Dict) -> Tuple[Dict, str]:
        """
        Prepare data for basic scatter plot (X vs Y)

        Args:
            formdata: Contains fileA, fileB, folder parameters

        Returns:
            Tuple of (scatter_data, message)
        """
        fileA = formdata.get('fileA', '')
        fileB = formdata.get('fileB', '')
        folder = formdata.get('folder', 'binned')

        if not fileA or not fileB:
            return {}, 'Please select both files for comparison'

        try:
            data_a = self.load_json_file(f'data/{folder}/{fileA}')
            data_b = self.load_json_file(f'data/{folder}/{fileB}')

            # Combine gene lists
            all_genes = set(data_a.keys()) | set(data_b.keys())

            scatter_data = []
            for gene in all_genes:
                scatter_data.append({
                    'gene': gene,
                    'x': float(data_a.get(gene, 0)),
                    'y': float(data_b.get(gene, 0))
                })

            return {
                'data': scatter_data,
                'fileA': fileA,
                'fileB': fileB,
                'total_genes': len(scatter_data)
            }, f'Prepared scatter plot with {len(scatter_data)} genes'

        except Exception as e:
            return {}, f'Error loading data: {str(e)}'

    def prepare_colored_scatter_data(self, formdata: Dict) -> Tuple[Dict, str]:
        """
        Prepare data for scatter plot colored by gene type or gene list

        Args:
            formdata: Contains files, folder, and coloring parameters

        Returns:
            Tuple of (scatter_data_with_colors, message)
        """
        fileA = formdata.get('fileA', '')
        fileB = formdata.get('fileB', '')
        folder = formdata.get('folder', 'binned')
        color_by = formdata.get('color_by', 'gene_type')
        gene_list = formdata.get('gene_list', '')

        if not fileA or not fileB:
            return {}, 'Please select both files for comparison'

        try:
            data_a = self.load_json_file(f'data/{folder}/{fileA}')
            data_b = self.load_json_file(f'data/{folder}/{fileB}')

            # Load gene type information
            gene_type_map = self.load_json_file('mappers/name2type.json')

            # Load gene list if specified
            highlight_genes = set()
            if gene_list:
                highlight_genes = set(self.load_json_file(f'data/genelist/{gene_list}'))

            # Combine gene lists
            all_genes = set(data_a.keys()) | set(data_b.keys())

            scatter_data = []
            for gene in all_genes:
                gene_type = gene_type_map.get(gene, 'unknown')
                in_list = gene in highlight_genes if highlight_genes else False

                scatter_data.append({
                    'gene': gene,
                    'x': float(data_a.get(gene, 0)),
                    'y': float(data_b.get(gene, 0)),
                    'type': gene_type,
                    'highlighted': in_list
                })

            return {
                'data': scatter_data,
                'fileA': fileA,
                'fileB': fileB,
                'color_by': color_by,
                'total_genes': len(scatter_data)
            }, f'Prepared colored scatter plot with {len(scatter_data)} genes'

        except Exception as e:
            return {}, f'Error loading data: {str(e)}'

    def prepare_gene_type_comparison(self, formdata: Dict) -> Tuple[Dict, str]:
        """
        Prepare data for gene type distribution comparison across datasets

        Args:
            formdata: Contains file list and parameters

        Returns:
            Tuple of (gene_type_data, message)
        """
        files = formdata.getlist('files')
        folder = formdata.get('folder', 'binned')

        if not files:
            return {}, 'Please select at least one file'

        try:
            gene_type_map = self.load_json_file('mappers/name2type.json')
            type_counts = {}

            for file in files:
                data = self.load_json_file(f'data/{folder}/{file}')
                file_type_counts = {}

                for gene, count in data.items():
                    gene_type = gene_type_map.get(gene, 'unknown')
                    file_type_counts[gene_type] = file_type_counts.get(gene_type, 0) + float(count)

                type_counts[file] = file_type_counts

            return {
                'type_counts': type_counts,
                'files': files
            }, f'Prepared gene type comparison for {len(files)} datasets'

        except Exception as e:
            return {}, f'Error loading data: {str(e)}'

    def prepare_heatmap_data(self, formdata: Dict) -> Tuple[Dict, str]:
        """
        Prepare data for expression heatmap

        Args:
            formdata: Contains file list, gene list, and parameters

        Returns:
            Tuple of (heatmap_data, message)
        """
        files = formdata.getlist('files')
        gene_list = formdata.get('gene_list', '')
        folder = formdata.get('folder', 'binned')

        if not files:
            return {}, 'Please select at least one file'

        try:
            # Load genes of interest
            if gene_list:
                genes_of_interest = set(self.load_json_file(f'data/genelist/{gene_list}'))
            else:
                # Use all genes from first file
                first_data = self.load_json_file(f'data/{folder}/{files[0]}')
                genes_of_interest = set(first_data.keys())

            # Prepare matrix
            expression_matrix = []
            for gene in genes_of_interest:
                gene_row = {'gene': gene, 'values': []}
                for file in files:
                    data = self.load_json_file(f'data/{folder}/{file}')
                    gene_row['values'].append(float(data.get(gene, 0)))
                expression_matrix.append(gene_row)

            return {
                'matrix': expression_matrix,
                'files': files,
                'genes': list(genes_of_interest),
                'total_genes': len(genes_of_interest)
            }, f'Prepared heatmap with {len(genes_of_interest)} genes x {len(files)} samples'

        except Exception as e:
            return {}, f'Error loading data: {str(e)}'

    def load_json_file(self, filepath: str) -> Dict:
        """
        Load and parse a JSON file

        Args:
            filepath: Path to JSON file

        Returns:
            Parsed JSON data

        Raises:
            FileNotFoundError: If file doesn't exist
            json.JSONDecodeError: If file isn't valid JSON
        """
        if not os.path.exists(filepath):
            raise FileNotFoundError(f'File not found: {filepath}')

        with open(filepath, 'r') as f:
            return json.load(f)

    def calculate_fold_change(self, value_a: float, value_b: float, pseudocount: float = 1.0) -> float:
        """
        Calculate fold change with pseudocount to avoid division by zero

        Args:
            value_a: Numerator value
            value_b: Denominator value
            pseudocount: Small value added to avoid division by zero

        Returns:
            Fold change (log2)
        """
        import math
        return math.log2((value_a + pseudocount) / (value_b + pseudocount))

    def normalize_data(self, data: Dict[str, float], method: str = 'total') -> Dict[str, float]:
        """
        Normalize expression data

        Args:
            data: Dictionary of gene -> expression value
            method: Normalization method ('total', 'median', 'quantile')

        Returns:
            Normalized data dictionary
        """
        if method == 'total':
            total = sum(data.values())
            if total > 0:
                return {gene: (value / total) * 1e6 for gene, value in data.items()}
            return data

        elif method == 'median':
            import statistics
            values = [v for v in data.values() if v > 0]
            if values:
                median_val = statistics.median(values)
                return {gene: (value / median_val) for gene, value in data.items()}
            return data

        else:
            # Default: return as-is
            return data