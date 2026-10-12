import React, { useEffect, useMemo, useRef } from 'react';
import * as d3 from 'd3';
import { MetricPlotData, MetricDataPoint } from '../types';

interface MetricsPlotProps {
  data: MetricPlotData | {};
  height: number;
  title: string;
  xAxisLabel: string;
  yAxisLabel: string;
  datasetId?: string;
  algorithm?: string;
}

const COLORS: Record<string, string> = {
  proportional: '#D76364',
  dense: '#5F97D2',
  fixed_spotDistance: '#965478',
  fixed_spotSize: '#8A8A88',
};

/**
 * Compute x-axis tick values matching the Python plotting logic.
 *
 * resolution_type: "spotsize" | "spotdistance"
 * simulations: the keys present in `data` (arrangement names)
 */
function computeXTicks(
  datasetId: string,
  algorithm: string,
  xAxisLabel: string,
  simulations: string[],
  metric: string,
): number[] | null {
  const isXe   = datasetId.startsWith('Xe');
  const isViHD = datasetId.startsWith('ViHD');
  if (!isXe && !isViHD) return null;

  const resType = xAxisLabel.toLowerCase().includes('distance') ? 'spotdistance' : 'spotsize';
  const isParnorm = metric === 'PARnorm';
  const simName = simulations.length === 1 ? simulations[0] : '';

  const arange = (start: number, stop: number, step: number): number[] => {
    const result: number[] = [];
    for (let v = start; v < stop; v += step) result.push(Math.round(v * 1e9) / 1e9);
    return result;
  };

  if (isXe) {
    if (isParnorm && (simName === 'proportional' || simName === 'dense')) {
      return resType === 'spotsize'
        ? arange(20, 151, 10)
        : arange(20, 281, 20);
    }
    if (isParnorm && simName === 'fixed_spotDistance') return arange(20, 101, 10);
    if (isParnorm && simName === 'fixed_spotSize')     return arange(55, 210, 15);
    // non-PARnorm
    if (resType === 'spotsize') return arange(20, 151, 10);
    return algorithm === 'TESLA' ? arange(55, 210, 15) : arange(20, 281, 20);
  }

  if (isViHD) {
    if (isParnorm && (simName === 'proportional' || simName === 'dense')) {
      return resType === 'spotsize'
        ? arange(16, 145, 16)
        : arange(16, 167, 32);
    }
    if (isParnorm && simName === 'fixed_spotDistance') return arange(16, 97, 8);
    if (isParnorm && simName === 'fixed_spotSize')     return arange(56, 217, 16);
    // non-PARnorm
    if (resType === 'spotsize') return arange(16, 145, 16);
    return algorithm === 'TESLA' ? arange(56, 217, 16) : arange(16, 273, 32);
  }

  return null;
}

const MetricsPlot: React.FC<MetricsPlotProps> = ({
  data,
  height,
  title,
  xAxisLabel,
  yAxisLabel,
  datasetId = '',
  algorithm = '',
}) => {
  const svgRef = useRef<SVGSVGElement>(null);

  const { totalWidth, w, margin, xMin, xMax } = useMemo(() => {
    const margin = { top: 40, right: 120, bottom: 60, left: 60 };

    const allPoints: MetricDataPoint[] = Object.values(data).flat();
    if (allPoints.length === 0) {
      return { totalWidth: 600, w: 420, margin, xMin: 0, xMax: 100 };
    }

    const xMin = d3.min(allPoints, d => d.axisValue) || 0;
    const xMax = d3.max(allPoints, d => d.axisValue) || 100;

    const pixelsPerUnit = 2;
    const chartAreaWidth = Math.max((xMax - xMin) * pixelsPerUnit, 400);
    const totalWidth = chartAreaWidth + margin.left + margin.right;

    return { totalWidth, w: chartAreaWidth, margin, xMin, xMax };
  }, [data]);

  useEffect(() => {
    if (!data || !svgRef.current) return;

    const allPoints: MetricDataPoint[] = Object.values(data).flat();
    if (allPoints.length === 0) return;

    const svg = d3.select(svgRef.current);
    svg.selectAll('*').remove();

    const h = height - margin.top - margin.bottom;
    const g = svg.append('g').attr('transform', `translate(${margin.left},${margin.top})`);

    // Y domain
    const yMin = d3.min(allPoints, d => d.metricValue) ?? 0;
    const yMax = d3.max(allPoints, d => d.metricValue) ?? 1;
    const yPad = (yMax - yMin) * 0.05 || 0.05;

    // X ticks — dataset-aware
    const simulations = Object.keys(data).filter(
      k => (data as Record<string, MetricDataPoint[]>)[k]?.length > 0,
    );
    const datasetTicks = computeXTicks(datasetId, algorithm, xAxisLabel, simulations, yAxisLabel);

    let tickValues: number[];
    let xDomainMin: number;
    let xDomainMax: number;

    if (datasetTicks && datasetTicks.length >= 2) {
      tickValues = datasetTicks;
      const step = tickValues[1] - tickValues[0];
      xDomainMin = tickValues[0] - step * 0.3;
      xDomainMax = tickValues[tickValues.length - 1] + step * 0.3;
    } else {
      // Fallback: derive ticks from data, step = 10
      const tickStart = Math.ceil(xMin / 10) * 10;
      const tickEnd   = Math.floor(xMax / 10) * 10;
      tickValues = d3.range(tickStart, tickEnd + 1, 10);
      const step = tickValues.length > 1 ? tickValues[1] - tickValues[0] : 10;
      xDomainMin = tickValues[0] - step * 0.3;
      xDomainMax = tickValues[tickValues.length - 1] + step * 0.3;
    }

    const xScale = d3.scaleLinear().domain([xDomainMin, xDomainMax]).range([0, w]);
    const yScale = d3.scaleLinear().domain([yMin - yPad, yMax + yPad]).range([h, 0]);

    // Axes — half-length ticks (3 px instead of default 6 px)
    const xAxis = d3.axisBottom(xScale)
      .tickValues(tickValues)
      .tickFormat(d => d.valueOf().toString())
      .tickSizeInner(3)
      .tickSizeOuter(0);

    const yAxis = d3.axisLeft(yScale)
      .ticks(5)
      .tickSizeInner(3)
      .tickSizeOuter(0);

    g.append('g')
      .attr('transform', `translate(0,${h})`)
      .call(xAxis)
      .selectAll<SVGTextElement, unknown>('text')
      .style('text-anchor', 'end')
      .attr('transform', 'rotate(-45)')
      .attr('dx', '-0.4em')
      .attr('dy', '0.6em');

    g.append('g').call(yAxis);

    // Title
    g.append('text')
      .attr('x', w / 2)
      .attr('y', -15)
      .attr('text-anchor', 'middle')
      .style('font-size', '14px')
      .style('font-weight', 'bold')
      .text(title);

    // X label
    g.append('text')
      .attr('x', w / 2)
      .attr('y', h + 55)
      .attr('text-anchor', 'middle')
      .style('font-size', '12px')
      .text(xAxisLabel);

    // Y label
    g.append('text')
      .attr('transform', 'rotate(-90)')
      .attr('y', -40)
      .attr('x', -h / 2)
      .attr('text-anchor', 'middle')
      .style('font-size', '12px')
      .text(yAxisLabel);

    const line = d3.line<MetricDataPoint>()
      .x(d => xScale(d.axisValue))
      .y(d => yScale(d.metricValue));

    Object.entries(data).forEach(([key, seriesData]) => {
      if (!seriesData || seriesData.length === 0) return;
      const color = COLORS[key] || '#000';

      // Line path
      g.append('path')
        .datum(seriesData)
        .attr('fill', 'none')
        .attr('stroke', color)
        .attr('stroke-width', 2)
        .attr('d', line);

      // Points
      seriesData.forEach((d: MetricDataPoint) => {
        g.append('circle')
          .attr('cx', xScale(d.axisValue))
          .attr('cy', yScale(d.metricValue))
          .attr('r', 3.5)
          .attr('fill', color);
      });
    });

    // Legend — only show series that have data
    const legend = svg.append('g')
      .attr('transform', `translate(${totalWidth - 110}, ${margin.top})`);

    let legendIndex = 0;
    Object.entries(data).forEach(([key, series]) => {
      if (!series || (series as MetricDataPoint[]).length === 0) return;
      const row = legend.append('g').attr('transform', `translate(0, ${legendIndex * 20})`);
      row.append('rect')
        .attr('width', 10)
        .attr('height', 10)
        .attr('fill', COLORS[key] || '#000');
      row.append('text')
        .attr('x', 15)
        .attr('y', 10)
        .style('font-size', '10px')
        .text(key);
      legendIndex++;
    });
  }, [data, height, title, xAxisLabel, yAxisLabel, totalWidth, w, margin, xMin, xMax, datasetId, algorithm]);

  return (
    <svg
      ref={svgRef}
      width={totalWidth}
      height={height}
      className="bg-white border border-gray-100 rounded shadow-sm"
    />
  );
};

export default MetricsPlot;
