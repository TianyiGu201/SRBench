import React, { useRef, useEffect, useMemo, useState } from 'react';
import * as d3 from 'd3';
import { SpatialPoint, GridRegion } from '../types';

interface SpatialViewerProps {
  imageUrl: string;
  width?: number;
  height?: number;
  points?: SpatialPoint[];
  regions?: GridRegion[];
  colorScaleType?: 'red-blue' | 'viridis';
  title?: string;
  onRegionClick?: (region: GridRegion) => void;
  selectedRegionId?: string | null;
  expressionDomain?: [number, number];
  showBackground?: boolean;
  shrinkRegions?: boolean;
}

const SpatialViewer: React.FC<SpatialViewerProps> = ({
  imageUrl,
  width = 400,
  height,
  points,
  regions,
  colorScaleType = 'viridis',
  title,
  onRegionClick,
  selectedRegionId,
  expressionDomain,
  showBackground = true,
  shrinkRegions = false,
}) => {
  const svgRef = useRef<SVGSVGElement>(null);
  const imgRef = useRef<HTMLImageElement>(null);

  const [imageNaturalSize, setImageNaturalSize] = useState({
    width: 1,
    height: 1,
  });

  const [displayHeight, setDisplayHeight] = useState<number>(height || 400);

  const handleImageLoad = () => {
    if (!imgRef.current) return;

    const naturalWidth = imgRef.current.naturalWidth;
    const naturalHeight = imgRef.current.naturalHeight;

    if (naturalWidth > 0 && naturalHeight > 0) {
      setImageNaturalSize({
        width: naturalWidth,
        height: naturalHeight,
      });

      if (!height) {
        const aspectRatio = naturalHeight / naturalWidth;
        setDisplayHeight(width * aspectRatio);
      }
    }
  };

  useEffect(() => {
    if (height) {
      setDisplayHeight(height);
    }
  }, [height, imageUrl]);

  /**
   * 这里直接使用 jpg 图像本身的 natural size 作为坐标 domain。
   * 前提：CSV 中的 region 坐标已经被批量转换到 jpg 坐标系。
   */
  const xScale = useMemo(() => {
    return d3
      .scaleLinear()
      .domain([0, imageNaturalSize.width])
      .range([0, width]);
  }, [imageNaturalSize.width, width]);

  const yScale = useMemo(() => {
    return d3
      .scaleLinear()
      .domain([0, imageNaturalSize.height])
      .range([0, displayHeight]);
  }, [imageNaturalSize.height, displayHeight]);

  const regionRenderConstraint = useMemo(() => {
    if (!regions || regions.length === 0 || !shrinkRegions) {
      return {
        safeWidthPx: null as number | null,
        safeHeightPx: null as number | null,
      };
    }

    const centersX = Array.from(
      new Set(regions.map((r) => xScale((r.xStart + r.xEnd) / 2)))
    ).sort((a, b) => a - b);

    const centersY = Array.from(
      new Set(regions.map((r) => yScale((r.yStart + r.yEnd) / 2)))
    ).sort((a, b) => a - b);

    let minDx = Infinity;
    let minDy = Infinity;

    for (let i = 1; i < centersX.length; i++) {
      const diff = centersX[i] - centersX[i - 1];
      if (diff > 0) minDx = Math.min(minDx, diff);
    }

    for (let i = 1; i < centersY.length; i++) {
      const diff = centersY[i] - centersY[i - 1];
      if (diff > 0) minDy = Math.min(minDy, diff);
    }

    const gapFactor = 0.92;

    return {
      safeWidthPx: Number.isFinite(minDx) ? minDx * gapFactor : null,
      safeHeightPx: Number.isFinite(minDy) ? minDy * gapFactor : null,
    };
  }, [regions, shrinkRegions, xScale, yScale]);

  useEffect(() => {
    if (!svgRef.current) return;

    const svg = d3.select(svgRef.current);
    svg.selectAll('*').remove();

    const drawTrueRegionBorder = (
      region: GridRegion,
      className: string,
      stroke: string,
      strokeWidth: number
    ) => {
      const x0 = xScale(region.xStart);
      const x1 = xScale(region.xEnd);
      const y0 = yScale(region.yStart);
      const y1 = yScale(region.yEnd);

      const x = Math.min(x0, x1);
      const y = Math.min(y0, y1);
      const w = Math.abs(x1 - x0);
      const h = Math.abs(y1 - y0);

      svg
        .append('rect')
        .attr('class', className)
        .attr('x', x)
        .attr('y', y)
        .attr('width', w)
        .attr('height', h)
        .attr('fill', 'none')
        .attr('stroke', stroke)
        .attr('stroke-width', strokeWidth)
        .attr('pointer-events', 'none');
    };

    if (regions && regions.length > 0) {
      const parNormMeans = regions
        .map((r) => r.parNormMean)
        .filter((v): v is number => v !== undefined && !Number.isNaN(v));

      const minPar = parNormMeans.length > 0 ? Math.min(...parNormMeans) : 0;
      const maxPar = parNormMeans.length > 0 ? Math.max(...parNormMeans) : 1;

      const colorDomain =
        minPar === maxPar ? [minPar - 1e-9, maxPar + 1e-9] : [minPar, maxPar];

      const colorScale = d3
        .scaleSequential<number>()
        .domain(
          colorScaleType === 'red-blue'
            ? [colorDomain[1], colorDomain[0]]
            : colorDomain
        )
        .interpolator(
          colorScaleType === 'red-blue' ? d3.interpolateRdBu : d3.interpolateViridis
        );

      const rects = svg
        .selectAll<SVGRectElement, GridRegion>('rect.region-cell')
        .data(regions)
        .enter()
        .append('rect')
        .attr('class', 'region-cell')
        .attr('fill', (d) => colorScale(d.parNormMean))
        .style('cursor', 'pointer');

      if (shrinkRegions) {
        rects
          .attr('x', (d) => {
            const cx = xScale((d.xStart + d.xEnd) / 2);
            const rawWidthPx = Math.abs(xScale(d.xEnd) - xScale(d.xStart));
            const safeWidthPx = regionRenderConstraint.safeWidthPx ?? rawWidthPx;
            const drawWidthPx = Math.min(rawWidthPx, safeWidthPx);
            return cx - drawWidthPx / 2;
          })
          .attr('y', (d) => {
            const cy = yScale((d.yStart + d.yEnd) / 2);
            const rawHeightPx = Math.abs(yScale(d.yEnd) - yScale(d.yStart));
            const safeHeightPx = regionRenderConstraint.safeHeightPx ?? rawHeightPx;
            const drawHeightPx = Math.min(rawHeightPx, safeHeightPx);
            return cy - drawHeightPx / 2;
          })
          .attr('width', (d) => {
            const rawWidthPx = Math.abs(xScale(d.xEnd) - xScale(d.xStart));
            const safeWidthPx = regionRenderConstraint.safeWidthPx ?? rawWidthPx;
            return Math.min(rawWidthPx, safeWidthPx);
          })
          .attr('height', (d) => {
            const rawHeightPx = Math.abs(yScale(d.yEnd) - yScale(d.yStart));
            const safeHeightPx = regionRenderConstraint.safeHeightPx ?? rawHeightPx;
            return Math.min(rawHeightPx, safeHeightPx);
          })
          .attr('opacity', 0.9)
          .attr('stroke', (d) => (d.id === selectedRegionId ? 'cyan' : 'none'))
          .attr('stroke-width', (d) => (d.id === selectedRegionId ? 2 : 0));

        rects
          .on('mouseover', function (_event, d) {
            svg.selectAll('.hover-border').remove();
            drawTrueRegionBorder(d, 'hover-border', '#FFFF00', 2);
          })
          .on('mouseout', function () {
            svg.selectAll('.hover-border').remove();
          });

        if (selectedRegionId) {
          const selected = regions.find((r) => r.id === selectedRegionId);
          if (selected) {
            drawTrueRegionBorder(selected, 'selected-border', 'cyan', 2);
          }
        }
      } else {
        rects
          .attr('x', (d) => Math.min(xScale(d.xStart), xScale(d.xEnd)))
          .attr('y', (d) => Math.min(yScale(d.yStart), yScale(d.yEnd)))
          .attr('width', (d) => Math.abs(xScale(d.xEnd) - xScale(d.xStart)))
          .attr('height', (d) => Math.abs(yScale(d.yEnd) - yScale(d.yStart)))
          .attr('opacity', 0.6)
          .attr('stroke', (d) => (d.id === selectedRegionId ? 'yellow' : 'white'))
          .attr('stroke-width', (d) => (d.id === selectedRegionId ? 2 : 0.5));
      }

      rects
        .on('click', (_event, d) => {
          if (onRegionClick) onRegionClick(d);
        })
        .append('title')
        .text(
          (d) =>
            `Region: ${d.id}\nminx: ${d.xStart}\nmaxx: ${d.xEnd}\nminy: ${d.yStart}\nmaxy: ${d.yEnd}\nPARnorm: ${d.parNormMean.toFixed(3)}`
        );
    }

    if (points && points.length > 0) {
      const visiblePoints = points.filter((p) => p.expression > 0.1);

      let minExpression: number;
      let maxExpression: number;

      if (expressionDomain) {
        [minExpression, maxExpression] = expressionDomain;
      } else {
        const expressionValues = visiblePoints.map((p) => p.expression);
        minExpression = d3.min(expressionValues) ?? 0;
        maxExpression = d3.max(expressionValues) ?? 1;
      }

      const safeDomain =
        minExpression === maxExpression
          ? [minExpression - 1e-9, maxExpression + 1e-9]
          : [minExpression, maxExpression];

      const colorScale = d3
        .scaleSequential(
          colorScaleType === 'red-blue' ? d3.interpolateRdBu : d3.interpolateViridis
        )
        .domain(safeDomain);

      svg
        .selectAll<SVGCircleElement, SpatialPoint>('circle.expression-point')
        .data(visiblePoints)
        .enter()
        .append('circle')
        .attr('class', 'expression-point')
        .attr('cx', (d) => xScale(d.x))
        .attr('cy', (d) => yScale(d.y))
        .attr('r', 2)
        .attr('fill', (d) => colorScale(d.expression))
        .attr('opacity', 0.8);
    }
  }, [
    points,
    regions,
    width,
    displayHeight,
    colorScaleType,
    selectedRegionId,
    onRegionClick,
    xScale,
    yScale,
    expressionDomain,
    shrinkRegions,
    regionRenderConstraint,
  ]);

  return (
    <div className="flex flex-col items-center">
      {title && <h4 className="text-sm font-semibold mb-2 text-gray-700">{title}</h4>}

      <div
        className="relative border border-gray-300 shadow-sm bg-gray-100 overflow-hidden"
        style={{ width, height: displayHeight }}
      >
        {showBackground && (
          <img
            ref={imgRef}
            src={imageUrl}
            alt="HE Background"
            className="absolute top-0 left-0"
            style={{
              width,
              height: displayHeight,
              display: 'block',
            }}
            onLoad={handleImageLoad}
          />
        )}

        <svg
          ref={svgRef}
          width={width}
          height={displayHeight}
          className="absolute top-0 left-0"
        />
      </div>
    </div>
  );
};

export default SpatialViewer;