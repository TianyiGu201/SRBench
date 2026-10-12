import React, {useEffect, useRef, useState} from 'react';
import colorMap from '../assets/rdbu-r.json';
import {decodeStaticResponse, getSpatialArrayUrl} from '../services/api';

type GeneArray = {rows: number; cols: number; values: Float64Array};
const arrayCache = new Map<string, Promise<GeneArray>>();

function loadNpy(url: string): Promise<GeneArray> {
  if (!arrayCache.has(url)) {
    arrayCache.set(url, fetch(url).then(async rawResponse => {
      const response = await decodeStaticResponse(rawResponse);
      const buffer = await response.arrayBuffer();
      const view = new DataView(buffer);
      if (view.getUint8(0) !== 0x93 || new TextDecoder().decode(buffer.slice(1, 6)) !== 'NUMPY') {
        throw new Error('Invalid NPY file');
      }
      const major = view.getUint8(6);
      const headerSize = major === 1 ? view.getUint16(8, true) : view.getUint32(8, true);
      const headerStart = major === 1 ? 10 : 12;
      const header = new TextDecoder().decode(buffer.slice(headerStart, headerStart + headerSize));
      const shape = header.match(/'shape':\s*\((\d+),\s*(\d+)\)/);
      if (!shape || !/'descr':\s*'<f8'/.test(header) || !/'fortran_order':\s*False/.test(header)) {
        throw new Error('Unsupported NPY format');
      }
      const rows = Number(shape[1]);
      const cols = Number(shape[2]);
      const offset = headerStart + headerSize;
      if (offset + rows * cols * 8 > buffer.byteLength) throw new Error('Truncated NPY file');
      const values = offset % 8 === 0
        ? new Float64Array(buffer, offset, rows * cols)
        : Float64Array.from({length: rows * cols}, (_, i) => view.getFloat64(offset + i * 8, true));
      return {rows, cols, values};
    }));
  }
  return arrayCache.get(url)!;
}

interface Props {
  dataset: string;
  gene: string;
  type: 'reference' | 'prediction';
  vmin?: number;
  vmax?: number;
}

const StaticSpatialImage: React.FC<Props> = ({dataset, gene, type, vmin, vmax}) => {
  const canvasRef = useRef<HTMLCanvasElement>(null);
  const [array, setArray] = useState<GeneArray | null>(null);
  const [error, setError] = useState('');
  const url = getSpatialArrayUrl(dataset, gene, type);

  useEffect(() => {
    let active = true;
    setArray(null);
    setError('');
    loadNpy(url).then(data => {if (active) setArray(data);}).catch(err => {if (active) setError(String(err));});
    return () => {active = false;};
  }, [url]);

  useEffect(() => {
    if (!array || !canvasRef.current) return;
    const {rows, cols, values} = array;
    const canvas = canvasRef.current;
    canvas.width = cols;
    canvas.height = rows;
    const context = canvas.getContext('2d');
    if (!context) return;
    let autoMin = Infinity;
    let autoMax = -Infinity;
    for (const value of values) {
      if (Number.isFinite(value)) {
        if (value < autoMin) autoMin = value;
        if (value > autoMax) autoMax = value;
      }
    }
    const low = Number.isFinite(vmin) ? vmin! : autoMin;
    const high = Number.isFinite(vmax) ? vmax! : autoMax;
    const range = high > low ? high - low : 1;
    const image = context.createImageData(cols, rows);
    for (let i = 0; i < values.length; i++) {
      const value = values[i];
      const color = Number.isFinite(value)
        ? colorMap[Math.max(0, Math.min(255, Math.round((value - low) / range * 255)))]
        : [255, 255, 255];
      const j = i * 4;
      image.data[j] = color[0];
      image.data[j + 1] = color[1];
      image.data[j + 2] = color[2];
      image.data[j + 3] = 255;
    }
    context.putImageData(image, 0, 0);
  }, [array, vmin, vmax]);

  if (error) return <div className="p-3 text-sm text-red-700">{error}</div>;
  if (!array) return <div className="p-3 text-sm text-gray-500">Loading expression image...</div>;
  return <canvas ref={canvasRef} role="img" aria-label={`${type} ${gene}`} className="w-full h-auto" />;
};

export default StaticSpatialImage;
