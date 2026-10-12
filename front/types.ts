export type Algorithm = 'iSTAR' | 'XFuse' | 'TESLA' | 'scstGCN' | 'SuperST';
export type SpotArrangement = 'dense' | 'proportional' | 'fixed_spotSize' | 'fixed_spotDistance';
export type MetricType = 'PCC' | 'SCC' | 'KLD' | 'JSD' | 'SSIM' | 'PSNR';

export interface DatasetMeta {
  id: string;
  name: string;
  species: string;
  technology: 'Xenium' | 'CosMx SMI' | 'Visium HD';
  tissue: string;
  description: string;
  heImageUrl: string; // Placeholder for /data_root/.../he.png
}

export interface MetricVisData {
  gene: string;
  value: number;
}

export interface GridRegion {
  id: string;
  xStart: number;
  xEnd: number;
  yStart: number;
  yEnd: number;
  parNormMean: number;
}

export interface GeneMetric {
  gene: string;
  PCC: number;
  SCC: number;
  KLD: number;
  JSD: number;
  SSIM: number;
  PSNR: number;
}

// Spatial Point for visualization
export interface SpatialData {
  gene: string;
  points: SpatialPoint[];
  shape: [number, number];
}

// Spatial Point for visualization
export interface SpatialPoint {
  x: number;
  y: number;
  expression: number;
}

// Data structures for the metrics plot
export interface MetricDataPoint {
  axisValue: number;  // Represents spot size or distance
  metricValue: number; // Represents the metric value (e.g., PCC)
  stdev: number;
}

export type MetricPlotSeries = 'proportional' | 'dense' | 'fixed_spotDistance' | 'fixed_spotSize';

export interface MetricPlotData {
  // Data for plots where x-axis is "spot size"
  proportional_spotSize?: MetricDataPoint[];
  dense_spotSize?: MetricDataPoint[];
  fixed_spotDistance_spotSize?: MetricDataPoint[];
  
  // Data for plots where x-axis is "distance between centers of spots"
  proportional_distance?: MetricDataPoint[];
  dense_distance?: MetricDataPoint[];
  fixed_spotSize_distance?: MetricDataPoint[];
}

// Data structure for replicate tables
export interface ReplicateData {
  rep1: GeneMetric[];
  rep2: GeneMetric[];
  rep3: GeneMetric[];
}
