import React, {useState} from 'react';
import {getSpatialImageUrl, isStaticDataMode} from '../services/api';
import StaticSpatialImage from './StaticSpatialImage';

interface GeneResultCardProps {
  gene: string;
  dataset: string;
  heUrl: string;
  colorMode: 'auto' | 'manual';
}

const GeneResultCard: React.FC<GeneResultCardProps> = ({
  gene,
  dataset,
  heUrl,
  colorMode,
}) => {
  const [localMin, setLocalMin] = useState<number>(0);
  const [localMax, setLocalMax] = useState<number>(1);

  const effectiveMin = colorMode === 'manual' ? localMin : undefined;
  const effectiveMax = colorMode === 'manual' ? localMax : undefined;

  const referenceImageUrl = getSpatialImageUrl(dataset, gene, 'reference', effectiveMin, effectiveMax);
  const predictionImageUrl = getSpatialImageUrl(dataset, gene, 'prediction', effectiveMin, effectiveMax);

  return (
    <div className="bg-white p-6 rounded-xl shadow-sm border border-gray-200">
      {/* Header & Controls */}
      <div className="flex flex-col xl:flex-row justify-between items-start xl:items-center mb-6 gap-6">
        <div className="flex items-center">
          <div className="bg-blue-100 text-blue-800 text-xs font-bold px-2 py-1 rounded uppercase mr-2">Gene</div>
          <h3 className="text-xl font-bold text-gray-900">{gene}</h3>
        </div>

        {/* Color Range Controls - only shown in manual mode */}
        {colorMode === 'manual' && (
          <div className="flex flex-col gap-1 w-64">
            <div className="flex justify-between text-xs text-gray-500 font-medium">
              <span>Min</span>
              <span>Max</span>
            </div>
            <div className="flex items-center gap-2">
              <input
                type="number"
                step="0.1"
                value={localMin}
                onChange={e => setLocalMin(parseFloat(e.target.value))}
                className="w-16 border border-gray-300 rounded px-1 py-0.5 text-xs text-center focus:ring-blue-500 focus:border-blue-500"
              />
              <div className="h-4 flex-1 rounded-full border border-gray-200 relative overflow-hidden">
                <div
                  className="absolute inset-0"
                  style={{ background: 'linear-gradient(to right, #2166ac, #f7f7f7, #b2182b)' }}
                />
              </div>
              <input
                type="number"
                step="0.1"
                value={localMax}
                onChange={e => setLocalMax(parseFloat(e.target.value))}
                className="w-16 border border-gray-300 rounded px-1 py-0.5 text-xs text-center focus:ring-blue-500 focus:border-blue-500"
              />
            </div>
          </div>
        )}
      </div>

      {/* Image Viewers */}
      <div className="overflow-x-auto overflow-y-hidden">
        <div className="w-fit min-w-full mx-auto">
          <div className="flex gap-8 w-fit">
            {/* HE Image */}
            <div className="flex flex-col items-center shrink-0">
              <h4 className="text-sm font-semibold mb-2 text-gray-700">HE Image</h4>
              <div
                className="relative border border-gray-300 shadow-sm bg-gray-100 overflow-hidden"
                style={{ width: 300 }}
              >
                <img
                  src={heUrl}
                  alt="HE Background"
                  className="w-full h-auto object-contain"
                />
              </div>
            </div>

            {/* Reference Image */}
            <div className="flex flex-col items-center shrink-0">
              <h4 className="text-sm font-semibold mb-2 text-gray-700">Reference Expression</h4>
              <div
                className="relative border border-gray-300 shadow-sm bg-gray-100 overflow-hidden"
                style={{ width: 300 }}
              >
                {isStaticDataMode
                  ? <StaticSpatialImage dataset={dataset} gene={gene} type="reference" vmin={effectiveMin} vmax={effectiveMax} />
                  : <img src={referenceImageUrl} alt={`Reference ${gene}`} className="w-full h-auto object-contain" key={referenceImageUrl} />}
              </div>
            </div>

            {/* Prediction Image */}
            <div className="flex flex-col items-center shrink-0">
              <h4 className="text-sm font-semibold mb-2 text-gray-700">Prediction Expression</h4>
              <div
                className="relative border border-gray-300 shadow-sm bg-gray-100 overflow-hidden"
                style={{ width: 300 }}
              >
                {isStaticDataMode
                  ? <StaticSpatialImage dataset={dataset} gene={gene} type="prediction" vmin={effectiveMin} vmax={effectiveMax} />
                  : <img src={predictionImageUrl} alt={`Prediction ${gene}`} className="w-full h-auto object-contain" key={predictionImageUrl} />}
              </div>
            </div>
          </div>
        </div>
      </div>
    </div>
  );
};

export default GeneResultCard;
