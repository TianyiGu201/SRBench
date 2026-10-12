import React, { useState, useEffect, useRef } from 'react';
import { Download, Search, ArrowUpDown } from 'lucide-react';
import {
  fetchDatasets,
  fetchGridRegions,
  fetchRegionData,
  downloadRegionCSV,
} from '../services/api';
import { DatasetMeta, Algorithm, GridRegion, GeneMetric } from '../types';
import SpatialViewer from '../components/SpatialViewer';

const OutOfSamplePrediction: React.FC = () => {
  const [datasets, setDatasets] = useState<DatasetMeta[]>([]);
  const [selectedDataset, setSelectedDataset] = useState<string>('');
  const [algorithm, setAlgorithm] = useState<Algorithm>('iSTAR');
  const [regions, setRegions] = useState<GridRegion[]>([]);
  const [selectedRegion, setSelectedRegion] = useState<GridRegion | null>(null);
  const [regionData, setRegionData] = useState<GeneMetric[]>([]);
  const [tableSearch, setTableSearch] = useState('');
  const [sortConfig, setSortConfig] = useState<{
    key: keyof GeneMetric;
    direction: 'asc' | 'desc';
  } | null>(null);

  const regionsRequestRef = useRef(0);
  const regionDataRequestRef = useRef(0);

  useEffect(() => {
    fetchDatasets().then((d) => {
      const filtered = d.filter((item) => item.technology === 'Xenium');
      setDatasets(filtered);
      if (filtered.length > 0) setSelectedDataset(filtered[0].id);
    });
  }, []);

  useEffect(() => {
    if (!selectedDataset) return;

    const reqId = ++regionsRequestRef.current;
    regionDataRequestRef.current += 1;

    setRegions([]);
    setSelectedRegion(null);
    setRegionData([]);
    setTableSearch('');
    setSortConfig(null);

    fetchGridRegions(selectedDataset, algorithm)
      .then((data) => {
        if (reqId !== regionsRequestRef.current) return;
        setRegions(data);
      })
      .catch((err) => {
        if (reqId !== regionsRequestRef.current) return;
        console.error('Failed to fetch regions:', err);
        setRegions([]);
      });
  }, [selectedDataset, algorithm]);

  useEffect(() => {
    if (!selectedDataset || !selectedRegion) {
      setRegionData([]);
      return;
    }

    const reqId = ++regionDataRequestRef.current;

    setRegionData([]);
    setTableSearch('');
    setSortConfig(null);

    fetchRegionData(selectedDataset, algorithm, selectedRegion.id)
      .then((data) => {
        if (reqId !== regionDataRequestRef.current) return;
        setRegionData(data);
      })
      .catch((err) => {
        if (reqId !== regionDataRequestRef.current) return;
        console.error('Failed to fetch region data:', err);
        setRegionData([]);
      });
  }, [selectedRegion, selectedDataset, algorithm]);

  const handleSort = (key: keyof GeneMetric) => {
    let direction: 'asc' | 'desc' = 'asc';
    if (sortConfig && sortConfig.key === key && sortConfig.direction === 'asc') {
      direction = 'desc';
    }
    setSortConfig({ key, direction });
  };

  const sortedData = React.useMemo(() => {
    if (!sortConfig) return regionData;
    return [...regionData].sort((a, b) => {
      const aVal = a[sortConfig.key];
      const bVal = b[sortConfig.key];
      if (aVal < bVal) return sortConfig.direction === 'asc' ? -1 : 1;
      if (aVal > bVal) return sortConfig.direction === 'asc' ? 1 : -1;
      return 0;
    });
  }, [regionData, sortConfig]);

  const filteredData = sortedData.filter((d) =>
    d.gene.toLowerCase().includes(tableSearch.toLowerCase())
  );

  const currentDataset = datasets.find((d) => d.id === selectedDataset);

  const parNormValues = regions
    .map((r) => r.parNormMean)
    .filter((v) => v !== undefined && !isNaN(v));

  const globalMin = parNormValues.length > 0 ? Math.min(...parNormValues) : 0;
  const globalMax = parNormValues.length > 0 ? Math.max(...parNormValues) : 1;

  const comboKey = `${selectedDataset}__${algorithm}`;

  return (
    <div className="flex h-full">
      {/* Left: Visualization */}
      <div className="w-1/2 min-h-0 flex flex-col border-r border-gray-200 bg-white">
        <div className="bg-white border-b border-gray-200 shadow-sm z-20 flex flex-col">
          {/* Top Bar */}
          <div className="px-8 py-5 border-b border-gray-100">
            <div>
              <h1 className="text-2xl font-bold text-gray-900 tracking-tight">
                FOV-based out-of-sample prediction
              </h1>
              <p className="text-sm text-gray-500 mt-1">
                Each square corresponds to one FOV used for training. Colors represent
                PARnorm derived from gene-level metrics (PCC, SCC, JSD, KLD). Hover to
                highlight the FOV region. Click to view detailed gene-level metrics.
                FOVs which failed to predict are excluded.
              </p>
            </div>
          </div>

          {/* Filters */}
          <div className="px-8 py-4 bg-gray-50/50 flex flex-wrap gap-6 items-end border-b border-gray-200">
            <div className="flex flex-col gap-1.5 flex-1 min-w-[150px]">
              <label className="text-[11px] font-bold text-gray-500 uppercase tracking-wider">
                Dataset
              </label>
              <select
                className="w-full rounded-lg border border-gray-300 bg-white px-3 py-2 text-sm shadow-sm focus:border-blue-500 focus:ring-1 focus:ring-blue-500"
                value={selectedDataset}
                onChange={(e) => {
                  setSelectedRegion(null);
                  setRegionData([]);
                  setTableSearch('');
                  setSortConfig(null);
                  regionDataRequestRef.current += 1;
                  setSelectedDataset(e.target.value);
                }}
              >
                {datasets.map((d) => (
                  <option key={d.id} value={d.id}>
                    {d.name}
                  </option>
                ))}
              </select>
            </div>

            <div className="flex flex-col gap-1.5 w-40">
              <label className="text-[11px] font-bold text-gray-500 uppercase tracking-wider">
                Method
              </label>
              <select
                className="w-full rounded-lg border border-gray-300 bg-white px-3 py-2 text-sm shadow-sm focus:border-blue-500 focus:ring-1 focus:ring-blue-500"
                value={algorithm}
                onChange={(e) => {
                  setSelectedRegion(null);
                  setRegionData([]);
                  setTableSearch('');
                  setSortConfig(null);
                  regionDataRequestRef.current += 1;
                  setAlgorithm(e.target.value as Algorithm);
                }}
              >
                {['XFuse', 'TESLA', 'iSTAR', 'scstGCN', 'SuperST_t', 'SuperST_t_g'].map(
                  (a) => (
                    <option key={a} value={a}>
                      {a}
                    </option>
                  )
                )}
              </select>
            </div>
          </div>
        </div>

        <div className="flex-1 min-h-0 overflow-auto bg-gray-50">
          {currentDataset && (
            <div className="w-fit mx-auto space-y-3 pb-4">
              <SpatialViewer
                imageUrl={currentDataset.heImageUrl}
                width={450}
                regions={regions}
                onRegionClick={(region) => {
                  setRegionData([]);
                  setTableSearch('');
                  setSortConfig(null);
                  setSelectedRegion(region);
                }}
                selectedRegionId={selectedRegion?.id}
                colorScaleType="red-blue"
                title="Click on a grid region to see metrics"
                shrinkRegions={true}
              />

              {/* Colorbar */}
              <div className="w-[450px] flex flex-col items-center gap-1">
                <div className="w-full flex items-center gap-2">
                  <span className="text-xs text-gray-500 w-12 text-right">
                    {globalMin.toFixed(2)}
                  </span>
                  <div
                    className="h-3 flex-1 rounded border border-gray-200"
                    style={{
                      background:
                        'linear-gradient(to right, #b2182b, #f7f7f7, #2166ac)',
                    }}
                  />
                  <span className="text-xs text-gray-500 w-12 text-left">
                    {globalMax.toFixed(2)}
                  </span>
                </div>
                <div className="w-full flex justify-between px-14">
                  <div className="flex flex-col items-start">
                    <span className="text-xs font-semibold text-red-700">
                      Lower PARnorm
                    </span>
                    <span className="text-xs text-gray-500">Better Prediction</span>
                  </div>
                  <div className="flex flex-col items-end">
                    <span className="text-xs font-semibold text-blue-700">
                      Higher PARnorm
                    </span>
                    <span className="text-xs text-gray-500">Worse Prediction</span>
                  </div>
                </div>
              </div>
            </div>
          )}
        </div>
      </div>

      {/* Right: Data Table */}
      <div key={comboKey} className="w-1/2 flex flex-col bg-white">
        <div className="p-4 border-b border-gray-200 flex justify-between items-center bg-gray-50">
          <h3 className="font-bold text-gray-800">
            {selectedRegion ? `Region: ${selectedRegion.id}` : 'Select a region'}
          </h3>
          {selectedRegion && (
            <button
              className="flex items-center px-3 py-1.5 bg-white border border-gray-300 rounded hover:bg-gray-50 text-sm"
              onClick={() =>
                downloadRegionCSV(selectedDataset, algorithm, selectedRegion.id)
              }
            >
              <Download size={14} className="mr-2" /> Download CSV
            </button>
          )}
        </div>

        <div className="p-2 border-b border-gray-200">
          <div className="relative">
            <Search className="absolute left-3 top-2.5 text-gray-400" size={16} />
            <input
              type="text"
              placeholder="Search gene..."
              className="w-full pl-9 pr-3 py-2 border border-gray-300 rounded text-sm focus:ring-blue-500 focus:border-blue-500"
              value={tableSearch}
              onChange={(e) => setTableSearch(e.target.value)}
            />
          </div>
        </div>

        <div className="flex-1 overflow-auto">
          <table className="min-w-full text-sm text-left">
            <thead className="text-xs text-gray-700 uppercase bg-gray-100 sticky top-0">
              <tr>
                <th
                  className="px-4 py-3 cursor-pointer hover:bg-gray-200"
                  onClick={() => handleSort('gene')}
                >
                  Gene
                </th>
                {['PCC', 'SCC', 'KLD', 'JSD', 'SSIM', 'PSNR'].map((m) => (
                  <th
                    key={m}
                    className="px-4 py-3 cursor-pointer hover:bg-gray-200 whitespace-nowrap"
                    onClick={() => handleSort(m as keyof GeneMetric)}
                  >
                    <div className="flex items-center">
                      {m}
                      {sortConfig?.key === m && (
                        <ArrowUpDown size={12} className="ml-1" />
                      )}
                    </div>
                  </th>
                ))}
              </tr>
            </thead>
            <tbody className="divide-y divide-gray-100">
              {filteredData.length > 0 ? (
                filteredData.map((row) => (
                  <tr key={row.gene} className="hover:bg-blue-50 transition-colors">
                    <td className="px-4 py-2 font-medium text-gray-900">{row.gene}</td>
                    <td className="px-4 py-2">{row.PCC?.toFixed(3)}</td>
                    <td className="px-4 py-2">{row.SCC?.toFixed(3)}</td>
                    <td className="px-4 py-2">{row.KLD?.toFixed(3)}</td>
                    <td className="px-4 py-2">{row.JSD?.toFixed(3)}</td>
                    <td className="px-4 py-2">{row.SSIM?.toFixed(3)}</td>
                    <td className="px-4 py-2">{row.PSNR?.toFixed(2)}</td>
                  </tr>
                ))
              ) : (
                <tr>
                  <td colSpan={7} className="text-center py-10 text-gray-500">
                    {selectedRegion
                      ? 'No matching genes found.'
                      : 'Click a region on the left to load data.'}
                  </td>
                </tr>
              )}
            </tbody>
          </table>
        </div>
      </div>
    </div>
  );
};

export default OutOfSamplePrediction;