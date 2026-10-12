import React, {useEffect, useRef, useState} from 'react';
import {Check, ChevronDown, Palette, Search} from 'lucide-react';
import {fetchDatasets, fetchMetricOptions, fetchSimilarGenes} from '../services/api';
import {DatasetMeta, MetricType} from '../types';
import GeneResultCard from '../components/GeneResultCard';

const MetricVisualization: React.FC = () => {
  // Data State
  const [datasets, setDatasets] = useState<DatasetMeta[]>([]);
  const [selectedDataset, setSelectedDataset] = useState<string>('');
  const [selectedMetric, setSelectedMetric] = useState<MetricType>('PCC');
  const [availableValues, setAvailableValues] = useState<string[]>([]);
  const [metricValue, setMetricValue] = useState<string>('');
  const [results, setResults] = useState<{gene: string, dataset: string, heUrl: string}[]>([]);
  const [loading, setLoading] = useState(false);

  // Color mode: 'auto' = per-gene auto scale; 'manual' = each card has its own min/max inputs
  const [colorMode, setColorMode] = useState<'auto' | 'manual'>('auto');

  // UI State
  const [isColorMenuOpen, setIsColorMenuOpen] = useState(false);
  const colorMenuRef = useRef<HTMLDivElement>(null);

  useEffect(() => {
    fetchDatasets().then(data => {
      setDatasets(data);
      if (data.length > 0) setSelectedDataset(data[0].id);
    });
  }, []);

  useEffect(() => {
    if (selectedDataset && selectedMetric) {
      fetchMetricOptions(selectedDataset, selectedMetric).then(opts => {
        setAvailableValues(opts);
          if (opts.length > 0) {
            const maxValue = Math.max(...opts.map(Number)).toString();
            setMetricValue(maxValue);
          }
      });
    }
  }, [selectedDataset, selectedMetric]);

  useEffect(() => {
    const handleClickOutside = (event: MouseEvent) => {
      if (colorMenuRef.current && !colorMenuRef.current.contains(event.target as Node)) {
        setIsColorMenuOpen(false);
      }
    };
    document.addEventListener('mousedown', handleClickOutside);
    return () => document.removeEventListener('mousedown', handleClickOutside);
  }, []);

  const handleVisualize = async () => {
    if (!selectedDataset || !metricValue) return;
    setLoading(true);
    setResults([]);
    try {
      const genes = await fetchSimilarGenes(selectedDataset, selectedMetric, metricValue);
      const dataset = datasets.find(d => d.id === selectedDataset);
      setResults(genes.map(gene => ({
        gene,
        dataset: selectedDataset,
        heUrl: dataset?.heImageUrl || '',
      })));
    } catch (error) {
      console.error('Error loading data:', error);
      alert('Error loading data. Please check the console for details.');
    } finally {
      setLoading(false);
    }
  };

  const SelectArrow = () => (
    <div className="pointer-events-none absolute inset-y-0 right-0 flex items-center px-2 text-gray-500">
      <ChevronDown size={14} />
    </div>
  );

  return (
    <div className="flex h-full flex-col bg-gray-50">
      {/* Header / Controls */}
      <div className="bg-white border-b border-gray-200 shadow-sm z-20 flex flex-col">

        {/* Top Bar */}
        <div className="px-8 py-5 border-b border-gray-100 flex flex-col md:flex-row md:items-center justify-between gap-4">
          <div>
            <h1 className="text-2xl font-bold text-gray-900 tracking-tight">Metric Visualization</h1>
            <p className="text-sm text-gray-500 mt-1">
              Select a dataset, metric, and target value to display the H&E image, ground truth expression, and predicted spatial expression for three genes whose metric values are closest to the selected target.
              Two color scale modes are available: Auto, which uses a shared color scale across all displayed genes, and Manual, which allows independent color ranges for each gene.
            </p>
          </div>
          <button
            onClick={handleVisualize}
            disabled={loading}
            className="flex items-center justify-center px-6 py-2.5 bg-blue-600 text-white text-sm font-semibold rounded-lg shadow-sm hover:bg-blue-700 focus:outline-none focus:ring-2 focus:ring-blue-500 focus:ring-offset-2 disabled:opacity-70 disabled:cursor-not-allowed transition-all min-w-[160px]"
          >
            {loading ? (
              <>
                <div className="animate-spin rounded-full h-4 w-4 border-2 border-white border-t-transparent mr-2"></div>
                Processing...
              </>
            ) : (
              <>
                <Search size={18} className="mr-2" />
                Visualize
              </>
            )}
          </button>
        </div>

        {/* Filters Toolbar */}
        <div className="px-8 py-4 bg-gray-50/50 flex flex-wrap gap-6 items-end border-b border-gray-200">

          {/* Dataset */}
          <div className="flex flex-col gap-1.5 flex-1 min-w-[200px]">
            <label className="text-[11px] font-bold text-gray-500 uppercase tracking-wider">Dataset</label>
            <div className="relative group">
              <select
                className="appearance-none pl-3 pr-8 py-2 w-full text-sm bg-white border border-gray-300 rounded-lg shadow-sm focus:border-blue-500 focus:ring-1 focus:ring-blue-500 cursor-pointer transition-shadow group-hover:border-gray-400"
                value={selectedDataset}
                onChange={(e) => setSelectedDataset(e.target.value)}
              >
                {datasets.map(d => <option key={d.id} value={d.id}>{d.name}</option>)}
              </select>
              <SelectArrow />
            </div>
          </div>

          {/* Metric */}
          <div className="flex flex-col gap-1.5 w-32">
            <label className="text-[11px] font-bold text-gray-500 uppercase tracking-wider">Metric</label>
            <div className="relative group">
              <select
                className="appearance-none pl-3 pr-8 py-2 w-full text-sm bg-white border border-gray-300 rounded-lg shadow-sm focus:border-blue-500 focus:ring-1 focus:ring-blue-500 cursor-pointer transition-shadow group-hover:border-gray-400"
                value={selectedMetric}
                onChange={(e) => setSelectedMetric(e.target.value as MetricType)}
              >
                {['PCC', 'SCC', 'KLD', 'JSD', 'SSIM', 'PSNR'].map(m => <option key={m} value={m}>{m}</option>)}
              </select>
              <SelectArrow />
            </div>
          </div>

          {/* Threshold */}
          <div className="flex flex-col gap-1.5 w-32">
            <label className="text-[11px] font-bold text-gray-500 uppercase tracking-wider">Value</label>
            <div className="relative group">
              <select
                className="appearance-none pl-3 pr-8 py-2 w-full text-sm bg-white border border-gray-300 rounded-lg shadow-sm focus:border-blue-500 focus:ring-1 focus:ring-blue-500 cursor-pointer transition-shadow group-hover:border-gray-400"
                value={metricValue}
                onChange={(e) => setMetricValue(e.target.value)}
              >
                {availableValues.map(v => <option key={v} value={v}>{v}</option>)}
              </select>
              <SelectArrow />
            </div>
          </div>

          {/* Color Scale Mode */}
          <div className="flex flex-col gap-1.5" ref={colorMenuRef}>
            <label className="text-[11px] font-bold text-gray-500 uppercase tracking-wider flex items-center gap-1">
              <Palette size={12} /> Color Scale
            </label>
            <div className="relative">
              <button
                onClick={() => setIsColorMenuOpen(!isColorMenuOpen)}
                className={`flex items-center gap-2 py-2 px-3 bg-white border border-gray-300 rounded-lg shadow-sm text-sm font-medium text-gray-700 hover:border-gray-400 focus:outline-none transition-all ${isColorMenuOpen ? 'border-blue-500 ring-1 ring-blue-500' : ''}`}
              >
                {colorMode === 'auto' ? 'Auto (Per Gene)' : 'Manual Range'}
                <ChevronDown size={14} className={`text-gray-400 transition-transform duration-200 ${isColorMenuOpen ? 'rotate-180' : ''}`} />
              </button>
              {isColorMenuOpen && (
                <div className="absolute top-full left-0 mt-2 w-44 bg-white rounded-lg shadow-xl border border-gray-100 py-1 z-50 ring-1 ring-black ring-opacity-5">
                  <button
                    onClick={() => { setColorMode('auto'); setIsColorMenuOpen(false); }}
                    className="w-full text-left px-4 py-3 text-sm hover:bg-gray-50 flex items-center justify-between transition-colors"
                  >
                    <span className={colorMode === 'auto' ? 'text-blue-600 font-medium' : 'text-gray-700'}>Auto (Per Gene)</span>
                    {colorMode === 'auto' && <Check size={16} className="text-blue-600" />}
                  </button>
                  <button
                    onClick={() => { setColorMode('manual'); setIsColorMenuOpen(false); }}
                    className="w-full text-left px-4 py-3 text-sm hover:bg-gray-50 flex items-center justify-between transition-colors"
                  >
                    <span className={colorMode === 'manual' ? 'text-blue-600 font-medium' : 'text-gray-700'}>Manual Range</span>
                    {colorMode === 'manual' && <Check size={16} className="text-blue-600" />}
                  </button>
                </div>
              )}
            </div>
          </div>

        </div>
      </div>

      {/* Results Area */}
      <div className="flex-1 overflow-auto p-8 bg-gray-50">
        {results.length === 0 && !loading && (
          <div className="flex flex-col items-center justify-center h-full text-gray-400 pb-20">
            <div className="bg-white p-6 rounded-full shadow-sm mb-4">
              <Search size={40} className="text-gray-300" />
            </div>
            <h3 className="text-lg font-medium text-gray-600">No results to display</h3>
            <p className="text-sm">Select a dataset and metric above, then click <b>Visualize</b>.</p>
          </div>
        )}
        <div className="space-y-12 max-w-7xl mx-auto min-w-0">
          {results.map((res, idx) => (
            <GeneResultCard
              key={res.gene + idx}
              gene={res.gene}
              dataset={res.dataset}
              heUrl={res.heUrl}
              colorMode={colorMode}
            />
          ))}
        </div>
      </div>
    </div>
  );
};

export default MetricVisualization;
