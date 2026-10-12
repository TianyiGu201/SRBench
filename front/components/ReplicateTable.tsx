import React from 'react';
import { GeneMetric, ReplicateData } from '../types';

interface ReplicateTableProps {
  data: ReplicateData | null;
}

const Table: React.FC<{ title: string; metrics: GeneMetric[] }> = ({ title, metrics }) => (
  <div className="flex-1">
    <h3 className="text-lg font-semibold text-gray-800 mb-2">{title}</h3>
    <div className="overflow-y-auto h-64 border border-gray-200 rounded-md">
      <table className="min-w-full divide-y divide-gray-200">
        <thead className="bg-gray-50">
          <tr>
            <th className="px-4 py-2 text-left text-xs font-medium text-gray-500 uppercase tracking-wider">Gene</th>
            <th className="px-4 py-2 text-left text-xs font-medium text-gray-500 uppercase tracking-wider">PCC</th>
            <th className="px-4 py-2 text-left text-xs font-medium text-gray-500 uppercase tracking-wider">SCC</th>
            <th className="px-4 py-2 text-left text-xs font-medium text-gray-500 uppercase tracking-wider">KLD</th>
            <th className="px-4 py-2 text-left text-xs font-medium text-gray-500 uppercase tracking-wider">JSD</th>
          </tr>
        </thead>
        <tbody className="bg-white divide-y divide-gray-200">
          {metrics.map((metric, index) => (
            <tr key={index}>
              <td className="px-4 py-2 whitespace-nowrap text-sm font-medium text-gray-900">{metric.gene}</td>
              <td className="px-4 py-2 whitespace-nowrap text-sm text-gray-500">{metric.PCC.toFixed(3)}</td>
              <td className="px-4 py-2 whitespace-nowrap text-sm text-gray-500">{metric.SCC.toFixed(3)}</td>
              <td className="px-4 py-2 whitespace-nowrap text-sm text-gray-500">{metric.KLD.toFixed(3)}</td>
              <td className="px-4 py-2 whitespace-nowrap text-sm text-gray-500">{metric.JSD.toFixed(3)}</td>
            </tr>
          ))}
        </tbody>
      </table>
    </div>
  </div>
);

const ReplicateTable: React.FC<ReplicateTableProps> = ({ data }) => {
  if (!data) {
    return <div className="text-center text-gray-500">Select parameters to load replicate data.</div>;
  }

  return (
    <div className="flex space-x-4">
      <Table title="Replicate 1" metrics={data.rep1} />
      <Table title="Replicate 2" metrics={data.rep2} />
      <Table title="Replicate 3" metrics={data.rep3} />
    </div>
  );
};

export default ReplicateTable;
