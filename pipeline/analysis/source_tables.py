"""Export the supplied gene-level Source Data without recomputing annotations."""
from __future__ import annotations
import csv
import json
from pathlib import Path


def export_workbook(workbook, output_dir, pipeline_gene_ids=False):
    from openpyxl import load_workbook
    source = Path(workbook)
    output = Path(output_dir)
    merged = output / 'merged'
    classified = output / 'classified'
    merged.mkdir(parents=True, exist_ok=True)
    classified.mkdir(parents=True, exist_ok=True)
    book = load_workbook(source, read_only=True, data_only=True)
    manifest = []
    gene_mapping = set()
    for sheet in book:
        if not sheet.title.endswith('_raw'):
            continue
        iterator = iter(sheet.values)
        metadata = {}
        header = None
        for row in iterator:
            if row and row[0] == 'gene' and 'PCC' in row:
                header = list(row)
                break
            if len(row) >= 2 and row[0] is not None:
                metadata[str(row[0])] = row[1]
        if header is None:
            continue
        rows = [dict(zip(header, row)) for row in iterator if row and row[0] is not None]
        dataset = metadata.get('Dataset name')
        method = metadata.get('Method')
        if not dataset or not method:
            raise ValueError(f'Missing dataset/method metadata in {sheet.title}')
        if len({r['gene'] for r in rows}) != len(rows):
            raise ValueError(f'Duplicate genes in {sheet.title}')
        for row in rows:
            row['dataset'] = dataset
            row['method'] = method
            if pipeline_gene_ids and str(dataset).startswith('cos_'):
                original = row['gene']
                row['source_gene'] = original
                row['gene'] = str(original).replace('/', '_')
                if row['gene'] != original:
                    gene_mapping.add((dataset, original, row['gene']))
        if len({row['gene'] for row in rows}) != len(rows):
            raise ValueError(f'Gene identifier mapping is not one-to-one in {sheet.title}')
        columns = [c for c in header if c is not None]
        for extra in ['dataset', 'method']:
            if extra not in columns:
                columns.append(extra)
        if pipeline_gene_ids and str(dataset).startswith('cos_'):
            columns.append('source_gene')
        stem = f'{dataset}_{method}'
        destinations = [merged / f'{stem}.csv']
        if {'hvg_status', 'morans_I', 'Spatial_Variability_Level', 'Gene_Expression_Level'}.issubset(columns):
            destinations.append(classified / f'{stem}_classified.csv')
        for dest in destinations:
            with dest.open('w', encoding='utf-8', newline='') as handle:
                writer = csv.DictWriter(handle, fieldnames=columns, extrasaction='ignore')
                writer.writeheader()
                writer.writerows(rows)
        manifest.append({'sheet': sheet.title, 'dataset': dataset, 'method': method,
                         'n_genes': len(rows), 'outputs': [p.relative_to(output).as_posix() for p in destinations]})
    book.close()
    if not manifest:
        raise ValueError('No *_raw gene-metric sheets found in the workbook')
    (output / 'export_manifest.json').write_text(json.dumps(manifest, indent=2) + '\n', encoding='utf-8')
    if pipeline_gene_ids:
        with (output / 'gene_identifier_mapping.csv').open('w', encoding='utf-8', newline='') as handle:
            writer = csv.writer(handle)
            writer.writerow(['dataset', 'source_gene', 'pipeline_gene'])
            writer.writerows(sorted(gene_mapping))
    print(f'Exported {len(manifest)} existing dataset-output tables to {output}')
    return manifest
