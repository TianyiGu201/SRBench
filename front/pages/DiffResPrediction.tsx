import React, { useState, useEffect, useMemo, useRef } from 'react';
import {
  fetchDatasets,
  fetchSpotArrangements,
  fetchResolutions,
  fetchMetricsPlotData,
  fetchReplicateData,
} from '../services/api';
import {
  DatasetMeta,
  Algorithm,
  SpotArrangement,
  MetricPlotData,
  ReplicateData,
  MetricType,
} from '../types';
import MetricsPlot from '../components/MetricsPlot';
import ReplicateTable from '../components/ReplicateTable';

const DiffResPrediction: React.FC = () => {
  const [datasets, setDatasets] = useState<DatasetMeta[]>([]);
  const [loading, setLoading] = useState(false);
  const [replicateLoading, setReplicateLoading] = useState(false);

  const [datasetId, setDatasetId] = useState('');
  const [algorithm, setAlgorithm] = useState<Algorithm>('iSTAR');
  const [metric, setMetric] = useState<MetricType>('PCC');
  const [plotData, setPlotData] = useState<MetricPlotData | null>(null);

  const [arrangements, setArrangements] = useState<string[]>([]);
  const [arrangement, setArrangement] = useState<SpotArrangement>('dense');
  const [resolutions, setResolutions] = useState<string[]>([]);
  const [resolution, setResolution] = useState('');
  const [replicateData, setReplicateData] = useState<ReplicateData | null>(null);

  // 用于防止旧请求覆盖新请求
  const plotRequestRef = useRef(0);
  const arrangementRequestRef = useRef(0);
  const resolutionRequestRef = useRef(0);
  const replicateRequestRef = useRef(0);

  // 初始化数据集
  useEffect(() => {
    let cancelled = false;

    fetchDatasets()
      .then((d) => {
        if (cancelled) return;

        const filtered = d.filter(
          (item) => item.technology === 'Xenium' || item.technology === 'Visium HD'
        );

        setDatasets(filtered);

        if (filtered.length > 0) {
          setDatasetId(filtered[0].id);
        } else {
          setDatasetId('');
        }
      })
      .catch((err) => {
        if (cancelled) return;
        console.error('Failed to fetch datasets:', err);
        setDatasets([]);
        setDatasetId('');
      });

    return () => {
      cancelled = true;
    };
  }, []);

  // dataset 或 algorithm 改变时，先清空下游状态，避免旧状态参与请求
  useEffect(() => {
    setPlotData(null);
    setArrangements([]);
    setResolutions([]);
    setResolution('');
    setReplicateData(null);

    if (algorithm === 'TESLA') {
      setArrangement('fixed_spotSize' as SpotArrangement);
    } else {
      setArrangement('' as SpotArrangement);
    }
  }, [datasetId, algorithm]);

  // arrangement 改变时，清空 resolution 和 replicateData
  useEffect(() => {
    setResolutions([]);
    setResolution('');
    setReplicateData(null);
  }, [datasetId, algorithm, arrangement]);

  // 请求 plot data
  useEffect(() => {
    if (!datasetId || !algorithm || !metric) return;

    const reqId = ++plotRequestRef.current;
    setLoading(true);
    setPlotData(null);

    fetchMetricsPlotData(datasetId, algorithm, metric)
      .then((data) => {
        if (reqId !== plotRequestRef.current) return;
        setPlotData(data);
      })
      .catch((err) => {
        if (reqId !== plotRequestRef.current) return;
        console.error('Failed to fetch plot data:', err);
        setPlotData(null);
      })
      .finally(() => {
        if (reqId !== plotRequestRef.current) return;
        setLoading(false);
      });
  }, [datasetId, algorithm, metric]);

  // 请求 arrangements
  useEffect(() => {
    if (!datasetId || !algorithm) return;

    const reqId = ++arrangementRequestRef.current;

    setArrangements([]);
    setResolutions([]);
    setResolution('');
    setReplicateData(null);

    if (algorithm === 'TESLA') {
      setArrangements(['fixed_spotSize']);
      setArrangement('fixed_spotSize' as SpotArrangement);
      return;
    }

    fetchSpotArrangements(datasetId, algorithm)
      .then((arr) => {
        if (reqId !== arrangementRequestRef.current) return;

        setArrangements(arr);
        if (arr.length > 0) {
          setArrangement(arr[0] as SpotArrangement);
        } else {
          setArrangement('' as SpotArrangement);
        }
      })
      .catch((err) => {
        if (reqId !== arrangementRequestRef.current) return;
        console.error('Failed to fetch spot arrangements:', err);
        setArrangements([]);
        setArrangement('' as SpotArrangement);
      });
  }, [datasetId, algorithm]);

  // 请求 resolutions
  useEffect(() => {
    if (!datasetId || !algorithm || !arrangement) return;

    const reqId = ++resolutionRequestRef.current;

    setResolutions([]);
    setResolution('');
    setReplicateData(null);

    fetchResolutions(datasetId, algorithm, arrangement)
      .then((res) => {
        if (reqId !== resolutionRequestRef.current) return;

        setResolutions(res);
        if (res.length > 0) {
          setResolution(res[0]);
        } else {
          setResolution('');
        }
      })
      .catch((err) => {
        if (reqId !== resolutionRequestRef.current) return;
        console.error('Failed to fetch resolutions:', err);
        setResolutions([]);
        setResolution('');
      });
  }, [datasetId, algorithm, arrangement]);

  // 请求 replicate data
  useEffect(() => {
    if (!datasetId || !algorithm || !arrangement || !resolution) return;

    const reqId = ++replicateRequestRef.current;

    setReplicateLoading(true);
    setReplicateData(null);

    fetchReplicateData(datasetId, algorithm, arrangement, resolution)
      .then((data) => {
        if (reqId !== replicateRequestRef.current) return;
        setReplicateData(data);
      })
      .catch((err) => {
        if (reqId !== replicateRequestRef.current) return;
        console.error('Failed to fetch replicate data:', err);
        setReplicateData(null);
      })
      .finally(() => {
        if (reqId !== replicateRequestRef.current) return;
        setReplicateLoading(false);
      });
  }, [datasetId, algorithm, arrangement, resolution]);

  const currentDataset = datasets.find((d) => d.id === datasetId);

  const { spotSizeData, distanceData } = useMemo(() => {
    if (!plotData) {
      return { spotSizeData: {}, distanceData: {} };
    }

    return {
      spotSizeData: {
        proportional: plotData.proportional_spotSize,
        dense: plotData.dense_spotSize,
        fixed_spotDistance: plotData.fixed_spotDistance_spotSize,
      },
      distanceData: {
        proportional: plotData.proportional_distance,
        dense: plotData.dense_distance,
        fixed_spotSize: plotData.fixed_spotSize_distance,
      },
    };
  }, [plotData]);

  return (
    <div className="flex flex-col h-full bg-gray-50 overflow-auto">
      <div className="bg-white border-b border-gray-200 p-6">
        <div className="mb-8">
          <h1 className="text-3xl font-bold text-gray-900 mb-2">
            Different Initial Resolution Prediction
          </h1>
          <p className="text-gray-600">
            Select a dataset, method, and evaluation metric to observe overall performance trends under various simulated resolutions. 
            The plots summarize the averaged results across simulations. 
            In the detailed section below, further select a specific simulation setting and resolution to examine gene-level prediction accuracy metrics across three replicates.
          </p>
        </div>

        <div className="grid grid-cols-1 md:grid-cols-2 gap-6 max-w-4xl">
          <div>
            <label className="block text-sm font-semibold text-gray-700 mb-2">
              DATASET
            </label>
            <select
              className="w-full rounded-md border-gray-300 shadow-sm focus:border-blue-500 focus:ring-blue-500 bg-gray-50 p-2.5 border"
              value={datasetId}
              onChange={(e) => setDatasetId(e.target.value)}
            >
              {datasets.map((d) => (
                <option key={d.id} value={d.id}>
                  {d.name}
                </option>
              ))}
            </select>
          </div>

          <div>
            <label className="block text-sm font-semibold text-gray-700 mb-2">
              METHOD
            </label>
            <select
              className="w-full rounded-md border-gray-300 shadow-sm focus:border-blue-500 focus:ring-blue-500 bg-gray-50 p-2.5 border"
              value={algorithm}
              onChange={(e) => setAlgorithm(e.target.value as Algorithm)}
            >
              {['XFuse', 'TESLA', 'iSTAR', 'scstGCN', 'SuperST_t', 'SuperST_t_g'].map((a) => (
                <option key={a} value={a}>
                  {a}
                </option>
              ))}
            </select>
          </div>
        </div>
      </div>

      <div className="p-8 space-y-10 max-w-7xl mx-auto w-full">
        <section>
          <div className="flex items-center justify-between mb-6">
            <h2 className="text-xl font-bold text-gray-800">Performance Visualization</h2>
            <div className="flex items-center">
              <span className="text-sm font-medium text-gray-600 mr-2">Metric:</span>
              <select
                className="rounded border-gray-300 text-sm p-1.5 bg-white border shadow-sm"
                value={metric}
                onChange={(e) => setMetric(e.target.value as MetricType)}
              >
                {['PCC', 'SCC', 'KLD', 'JSD'].map((m) => (
                  <option key={m} value={m}>
                    {m}
                  </option>
                ))}
              </select>
            </div>
          </div>

          <div className="flex flex-col gap-8 items-start">
            <div className="flex flex-col items-center">
              <h4 className="text-sm font-semibold mb-2 text-gray-700">HE Image</h4>
              <div
                className="relative border border-gray-300 shadow-sm bg-gray-100 overflow-hidden"
                style={{ width: 300 }}
              >
                <img
                  src={currentDataset?.heImageUrl}
                  alt="HE Background"
                  className="w-full h-auto object-contain"
                />
              </div>
            </div>

            <div className="w-full overflow-x-auto pb-4">
              {loading ? (
                <div className="h-64 flex items-center justify-center text-gray-400">
                  Loading charts...
                </div>
              ) : !plotData ? (
                <div className="h-64 flex items-center justify-center text-gray-400">
                  No chart data available
                </div>
              ) : (
                <div className="flex gap-8 min-w-max items-start">
                  {algorithm !== 'TESLA' && (
                    <MetricsPlot
                      data={spotSizeData}
                      height={350}
                      title="Initial Spatial Resolution of Different Spot Sizes"
                      xAxisLabel="Spot Size (μm)"
                      yAxisLabel={metric}
                      datasetId={datasetId}
                      algorithm={algorithm}
                    />
                  )}
                  <MetricsPlot
                    data={distanceData}
                    height={350}
                    title="Initial Spatial Resolution of Different Spot Distances"
                    xAxisLabel="Distance between Centers of Spots (μm)"
                    yAxisLabel={metric}
                    datasetId={datasetId}
                    algorithm={algorithm}
                  />
                </div>
              )}
            </div>
          </div>
        </section>

        <hr className="border-gray-200" />

        <section>
          <h2 className="text-xl font-bold text-gray-800 mb-4">Detailed Replicate Results</h2>

          <div className="flex items-end gap-6 mb-6">
            <div>
              <label className="block text-xs font-semibold text-gray-500 uppercase mb-1">
                Spot Layout
              </label>
              <select
                className="w-48 rounded border-gray-300 bg-white p-2 border text-sm disabled:bg-gray-100 disabled:text-gray-500"
                value={arrangement}
                onChange={(e) => setArrangement(e.target.value as SpotArrangement)}
                disabled={algorithm === 'TESLA'}
              >
                {arrangements.map((a) => (
                  <option key={a} value={a}>
                    {a}
                  </option>
                ))}
              </select>
            </div>

            <div>
              <label className="block text-xs font-semibold text-gray-500 uppercase mb-1">
                Initial Resolution
              </label>
              <select
                className="w-48 rounded border-gray-300 bg-white p-2 border text-sm"
                value={resolution}
                onChange={(e) => setResolution(e.target.value)}
              >
                {resolutions.map((r) => (
                  <option key={r} value={r}>
                    {r}
                  </option>
                ))}
              </select>
            </div>
          </div>

          {replicateLoading ? (
            <div className="text-center text-gray-500">Loading...</div>
          ) : !replicateData ? (
            <div className="text-center text-gray-400">No replicate data available</div>
          ) : (
            <ReplicateTable data={replicateData} />
          )}
        </section>
      </div>
    </div>
  );
};

export default DiffResPrediction;