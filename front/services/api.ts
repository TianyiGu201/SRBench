import {DatasetMeta, GridRegion, GeneMetric, SpatialPoint, SpatialData} from '../types';
import * as d3 from 'd3';

const API_BASE_URL = import.meta.env.VITE_API_BASE_URL || '';
const DATA_BASE_URL = (import.meta.env.VITE_DATA_BASE_URL || '').replace(/\/$/, '');
export const isStaticDataMode = Boolean(DATA_BASE_URL);
export const isCompressedDataMode = isStaticDataMode && import.meta.env.VITE_DATA_COMPRESSION === 'gzip';
const dataUrl = (path: string) => `${DATA_BASE_URL}${path}`;
const staticFileUrl = (path: string) => dataUrl(`${path}${isCompressedDataMode ? '.gz' : ''}`);

export async function decodeStaticResponse(response: Response): Promise<Response> {
  if (!response.ok) throw new Error(`Static data request failed: ${response.status} ${response.url}`);
  if (!isCompressedDataMode) return response;
  if (!response.body) throw new Error(`Empty compressed response: ${response.url}`);
  if (typeof DecompressionStream === 'undefined') {
    throw new Error('This browser does not support gzip decompression');
  }
  return new Response(response.body.pipeThrough(new DecompressionStream('gzip')));
}

async function fetchStaticJson<T>(path: string): Promise<T> {
  return (await decodeStaticResponse(await fetch(staticFileUrl(path)))).json() as Promise<T>;
}

export const getSpatialArrayUrl = (datasetId: string, gene: string, type: 'reference' | 'prediction') =>
  staticFileUrl(`/data_root/${encodeURIComponent(datasetId)}/metric_visualization/${type}/${encodeURIComponent(gene)}.npy`);

// Helper function for API calls
async function fetchAPI<T>(endpoint: string): Promise<T> {
  const response = await fetch(`${API_BASE_URL}${endpoint}`);
  if (!response.ok) {
    throw new Error(`API error: ${response.statusText}`);
  }
  return response.json();
}

export const fetchDatasets = async (): Promise<DatasetMeta[]> => {
  if (isStaticDataMode) {
    const datasets = await fetchStaticJson<DatasetMeta[]>('/site/datasets.json');
    return datasets.map(d => ({...d, heImageUrl: dataUrl(d.heImageUrl)}));
  }
  const datasets = await fetchAPI<DatasetMeta[]>('/api/datasets');
  // Convert relative image URLs to absolute URLs
  return datasets.map(d => ({
    ...d,
    heImageUrl: `${API_BASE_URL}${d.heImageUrl}`
  }));
};

export const fetchDatasetMetadata = async (datasetId: string) => {
  if (isStaticDataMode) {
    return fetchStaticJson(`/data_root/${encodeURIComponent(datasetId)}/dataset_metaINFO.json`);
  }
  return fetchAPI<any>(`/api/datasets/${datasetId}/metadata`);
};

export const fetchMetricOptions = async (datasetId: string, metric: string) => {
  if (isStaticDataMode) {
    const data = await fetchStaticJson<Record<string, Record<string, string[]>>>(`/data_root/${encodeURIComponent(datasetId)}/metric_vis.json`);
    return Object.keys(data[metric] || {});
  }
  const data = await fetchAPI<Record<string, Record<string, string[]>>>(`/api/metric-vis/${datasetId}/options`);
  // Return the keys for the specific metric
  if (data[metric]) {
    return Object.keys(data[metric]);
  }
  return [];
};

export const fetchSimilarGenes = async (datasetId: string, metric: string, value: string) => {
  if (isStaticDataMode) {
    const data = await fetchStaticJson<Record<string, Record<string, string[]>>>(`/data_root/${encodeURIComponent(datasetId)}/metric_vis.json`);
    const genes = data[metric]?.[value];
    if (!genes) throw new Error(`Metric value ${metric}/${value} is unavailable`);
    return genes as string[];
  }
  const data = await fetchAPI<{ genes: string[] }>(`/api/metric-vis/${datasetId}/${metric}/${value}/genes`);
  return data.genes;
};

export const fetchGridRegions = async (datasetId: string, algorithm: string): Promise<GridRegion[]> => {
  if (isStaticDataMode) {
    const data = await fetchStaticJson<{regions: GridRegion[]}>(`/site/oos/${encodeURIComponent(datasetId)}/${encodeURIComponent(algorithm)}/regions.json`);
    return data.regions;
  }
  const data = await fetchAPI<{ regions: GridRegion[] }>(`/api/out-of-sample/${datasetId}/${algorithm}/regions`);
  return data.regions;
};

export const fetchRegionData = async (datasetId: string, algorithm: string, regionId: string): Promise<GeneMetric[]> => {
  if (isStaticDataMode) {
    const data = await fetchStaticJson<{data: GeneMetric[]}>(`/site/oos/${encodeURIComponent(datasetId)}/${encodeURIComponent(algorithm)}/region/${encodeURIComponent(regionId)}.json`);
    return data.data;
  }
  const data = await fetchAPI<{
    region_id: string;
    data: GeneMetric[]
  }>(`/api/out-of-sample/${datasetId}/${algorithm}/region/${regionId}`);
  return data.data;
};

export const computeGlobalExpressionDomain = async (datasetId: string, gene: string) => {
  const data = await fetchAPI<[number, number]>(`/api/spatial-data/${datasetId}/compute-global-expression-domain/${gene}`);
  return data;
};

export const fetchSpatialData = async (datasetId: string, gene: string, type: 'reference' | 'prediction' = 'reference'): Promise<SpatialData> => {
  const endpoint = type === 'reference'
    ? `/api/spatial-data/${datasetId}/reference/${gene}`
    : `/api/spatial-data/${datasetId}/prediction/${gene}`;

  const data = await fetchAPI<SpatialData>(endpoint);
  return data;
};

// Get spatial image URL from backend (for matplotlib-generated images)
export const getSpatialImageUrl = (
  datasetId: string,
  gene: string,
  type: 'reference' | 'prediction' = 'reference',
  vmin?: number,
  vmax?: number
): string => {
  const params = new URLSearchParams();
  if (vmin !== undefined) params.append('vmin', vmin.toString());
  if (vmax !== undefined) params.append('vmax', vmax.toString());

  const queryString = params.toString();
  const endpoint = type === 'reference'
    ? `/api/spatial-image/${datasetId}/reference/${gene}`
    : `/api/spatial-image/${datasetId}/prediction/${gene}`;

  return `${API_BASE_URL}${endpoint}${queryString ? '?' + queryString : ''}`;
};

// New functions for different resolution prediction
export const fetchSpotArrangements = async (datasetId: string, algorithm: string): Promise<string[]> => {
  if (isStaticDataMode) {
    const manifest = await fetchStaticJson<{arrangements: string[]}>(`/site/diff/${encodeURIComponent(datasetId)}/${encodeURIComponent(algorithm)}/manifest.json`);
    return manifest.arrangements;
  }
  const data = await fetchAPI<{ arrangements: string[] }>(`/api/diff-res/${datasetId}/${algorithm}/arrangements`);
  return data.arrangements;
};

export const fetchResolutions = async (datasetId: string, algorithm: string, arrangement: string): Promise<string[]> => {
  if (isStaticDataMode) {
    const manifest = await fetchStaticJson<{resolutions: Record<string, string[]>}>(`/site/diff/${encodeURIComponent(datasetId)}/${encodeURIComponent(algorithm)}/manifest.json`);
    return manifest.resolutions[arrangement] || [];
  }
  const data = await fetchAPI<{
    resolutions: string[]
  }>(`/api/diff-res/${datasetId}/${algorithm}/${arrangement}/resolutions`);
  return data.resolutions;
};

export const fetchMetricsPlotData = async (datasetId: string, algorithm: string, metric: string) => {
  if (isStaticDataMode) {
    return fetchStaticJson(`/site/diff/${encodeURIComponent(datasetId)}/${encodeURIComponent(algorithm)}/plots/${encodeURIComponent(metric)}.json`);
  }
  return fetchAPI<any>(`/api/diff-res/${datasetId}/${algorithm}/${metric}/plot-data`);
};

export const fetchReplicateData = async (
  datasetId: string,
  algorithm: string,
  arrangement: string,
  resolution: string
) => {
  if (isStaticDataMode) {
    const result: Record<string, GeneMetric[]> = {};
    await Promise.all(['rep1', 'rep2', 'rep3'].map(async rep => {
      const path = `/data_root/${encodeURIComponent(datasetId)}/prediction/${encodeURIComponent(algorithm)}/raw_resolution/${rep}/${encodeURIComponent(arrangement)}/${encodeURIComponent(resolution)}.csv`;
      const response = await fetch(staticFileUrl(path));
      result[rep] = response.ok
        ? d3.csvParse(await (await decodeStaticResponse(response)).text(), d3.autoType) as unknown as GeneMetric[]
        : [];
    }));
    return result;
  }
  return fetchAPI<any>(`/api/diff-res/${datasetId}/${algorithm}/${arrangement}/${resolution}/replicates`);
};

// Helper to download files
export const downloadRegionCSV = (datasetId: string, algorithm: string, regionId: string) => {
  if (isStaticDataMode) {
    void fetchRegionData(datasetId, algorithm, regionId).then(data => {
      const csv = d3.csvFormat(data as any);
      const url = URL.createObjectURL(new Blob([csv], {type: 'text/csv'}));
      const link = document.createElement('a');
      link.href = url;
      link.download = `${datasetId}_${algorithm}_${regionId}.csv`;
      link.click();
      setTimeout(() => URL.revokeObjectURL(url), 1000);
    });
    return;
  }
  window.open(`${API_BASE_URL}/api/out-of-sample/${datasetId}/${algorithm}/region/${regionId}/download`, '_blank');
};

export const downloadHEImage = (datasetId: string) => {
  if (isStaticDataMode) {
    window.open(dataUrl(`/data_root/${encodeURIComponent(datasetId)}/he-raw.tiff`), '_blank');
    return;
  }
  window.open(`${API_BASE_URL}/api/datasets/${datasetId}/he-image`, '_blank');
};
