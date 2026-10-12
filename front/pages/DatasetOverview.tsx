import React, { useEffect, useState } from 'react';
import { Download, Database, Dna, Microscope, X } from 'lucide-react';
import { fetchDatasets, fetchDatasetMetadata, downloadHEImage } from '../services/api';
import { DatasetMeta } from '../types';

const DatasetOverview: React.FC = () => {
  const [datasets, setDatasets] = useState<DatasetMeta[]>([]);
  const [loading, setLoading] = useState(true);
  const [filterTech, setFilterTech] = useState<string>('All');
  const [selectedImage, setSelectedImage] = useState<string | null>(null);

  const [selectedMetadata, setSelectedMetadata] = useState<any | null>(null);
  const [metadataLoading, setMetadataLoading] = useState(false);

  useEffect(() => {
    fetchDatasets().then(data => {
      setDatasets(data);
      setLoading(false);
    });
  }, []);

  const filteredDatasets = filterTech === 'All'
    ? datasets
    : datasets.filter(d => d.technology === filterTech);

  const technologies = ['All', 'Xenium', 'Visium HD', 'CosMx SMI'];

  const openMetadata = async (datasetId: string) => {
    setMetadataLoading(true);
    try {
      const data = await fetchDatasetMetadata(datasetId);
      setSelectedMetadata(data);
    } catch (e) {
      console.error(e);
    } finally {
      setMetadataLoading(false);
    }
  };

  if (loading) return <div className="p-10 text-center">Loading datasets...</div>;

  return (
    <div className="p-8 max-w-7xl mx-auto">

      {/* HE Image Modal */}
      {selectedImage && (
        <div
          className="fixed inset-0 z-50 flex items-center justify-center bg-black/95 p-4 backdrop-blur-sm"
          onClick={() => setSelectedImage(null)}
        >
          <button
            className="absolute top-6 right-6 text-white/80 hover:text-white transition-colors"
            onClick={() => setSelectedImage(null)}
          >
            <X size={32} />
          </button>
          <img
            src={selectedImage}
            alt="Full Resolution HE"
            className="max-w-full max-h-full object-contain rounded shadow-2xl"
            onClick={(e) => e.stopPropagation()}
          />
        </div>
      )}

      {/* Metadata Modal */}
      {selectedMetadata && (
        <div
          className="fixed inset-0 z-50 flex items-center justify-center bg-black/80 p-6 backdrop-blur-sm"
          onClick={() => setSelectedMetadata(null)}
        >
          <div
            className="bg-white rounded-xl shadow-xl max-w-2xl w-full max-h-[80vh] overflow-auto p-6 relative"
            onClick={(e) => e.stopPropagation()}
          >
            <button
              className="absolute top-4 right-4 text-gray-500 hover:text-gray-800"
              onClick={() => setSelectedMetadata(null)}
            >
              <X size={24} />
            </button>

            <h3 className="text-xl font-bold mb-6">Dataset Metadata</h3>

            {metadataLoading ? (
              <div className="text-gray-500">Loading metadata...</div>
            ) : (
              <div className="grid grid-cols-2 gap-x-6 gap-y-3 text-sm">
                {Object.entries(selectedMetadata).map(([key, value]) => (
                  <React.Fragment key={key}>
                    <div className="font-medium text-gray-600">{key}</div>

                    <div className="text-gray-900 break-all">
                      {key === 'data_source_url' ? (
                        <a
                          href={String(value)}
                          target="_blank"
                          rel="noopener noreferrer"
                          className="text-blue-600 hover:underline"
                        >
                          {String(value)}
                        </a>
                      ) : (
                        String(value)
                      )}
                    </div>

                  </React.Fragment>
                ))}
              </div>
            )}
          </div>
        </div>
      )}

      <div className="mb-8">
        <h2 className="text-3xl font-bold text-gray-900 mb-2">Dataset Overview</h2>
        <p className="text-gray-600">Explore available spatial transcriptomics datasets and metadata.</p>
      </div>

      {/* Filter Tabs */}
      <div className="flex space-x-2 mb-8 border-b border-gray-200 pb-1">
        {technologies.map(tech => (
          <button
            key={tech}
            onClick={() => setFilterTech(tech)}
            className={`px-4 py-2 text-sm font-medium rounded-t-lg transition-colors ${
              filterTech === tech
                ? 'bg-blue-600 text-white shadow-sm'
                : 'bg-white text-gray-600 hover:bg-gray-100 hover:text-gray-900'
            }`}
          >
            {tech}
          </button>
        ))}
      </div>

      {/* Dataset Grid */}
      <div className="grid grid-cols-1 md:grid-cols-2 lg:grid-cols-3 gap-6">
        {filteredDatasets.map(dataset => (
          <div
            key={dataset.id}
            className="bg-white rounded-xl shadow-sm border border-gray-200 overflow-hidden hover:shadow-md transition-shadow group flex flex-col h-full"
          >
            {/* Header Image */}
            <div
              className="relative h-56 bg-gray-900 overflow-hidden cursor-pointer"
              onClick={() => setSelectedImage(dataset.heImageUrl)}
            >
              <img
                src={dataset.heImageUrl}
                alt={dataset.name}
                className="w-full h-full object-contain transition-transform duration-500 group-hover:scale-105 opacity-90 group-hover:opacity-100"
              />

              <div className="absolute inset-0 bg-black/40 opacity-0 group-hover:opacity-100 transition-opacity flex items-center justify-center">
                <span className="text-white text-sm font-medium bg-black/50 px-4 py-2 rounded-full border border-white/25 backdrop-blur-md flex items-center">
                  <Microscope size={16} className="mr-2" />
                  Click to view full HE image
                </span>
              </div>
            </div>

            {/* Content */}
            <div className="p-5 flex-1 flex flex-col">
              <div className="flex justify-between items-start mb-2">
                <h3 className="text-lg font-bold text-gray-900 line-clamp-1" title={dataset.name}>
                  {dataset.name}
                </h3>

                <span
                  className={`px-2 py-1 rounded text-xs font-semibold ${
                    dataset.technology === 'Xenium'
                      ? 'bg-purple-100 text-purple-700'
                      : dataset.technology === 'CosMx SMI'
                      ? 'bg-cyan-100 text-cyan-700'
                      : 'bg-orange-100 text-orange-700'
                  }`}
                >
                  {dataset.technology}
                </span>
              </div>

              <div className="space-y-2 mb-4 flex-1">
                <div className="flex items-center text-sm text-gray-600">
                  <Dna size={16} className="mr-2" />
                  <span>{dataset.species}</span>
                </div>

                <div className="flex items-center text-sm text-gray-600">
                  <Microscope size={16} className="mr-2" />
                  <span>{dataset.tissue}</span>
                </div>

                <p className="text-xs text-gray-500 line-clamp-3 mt-2">
                  {dataset.description}
                </p>
              </div>

              {/* Actions */}
              <div className="flex gap-2 mt-4 pt-4 border-t border-gray-100">
                <button
                  className="flex-1 flex items-center justify-center px-4 py-2 bg-gray-50 text-gray-700 text-sm font-medium rounded-lg hover:bg-gray-100 border border-gray-200 transition-colors"
                  onClick={(e) => {
                    e.stopPropagation();
                    downloadHEImage(dataset.id);
                  }}
                >
                  <Download size={16} className="mr-2" />
                  HE Image
                </button>

                <button
                  className="flex-1 flex items-center justify-center px-4 py-2 bg-gray-50 text-gray-700 text-sm font-medium rounded-lg hover:bg-gray-100 border border-gray-200 transition-colors"
                  onClick={(e) => {
                    e.stopPropagation();
                    openMetadata(dataset.id);
                  }}
                >
                  <Database size={16} className="mr-2" />
                  Metadata
                </button>
              </div>
            </div>
          </div>
        ))}
      </div>
    </div>
  );
};

export default DatasetOverview;
