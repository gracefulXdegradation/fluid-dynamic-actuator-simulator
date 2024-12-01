import React, { useEffect, useRef, useState } from "react";
import * as d3 from "d3";

const colors = ["steelblue", "orange", "green", "#ff6361"];

const height = (w: number) => w * 2 / 3 ;
const initialWidth = 600;

const LineGraph = ({ timestamps, values, graphNames }) => {
  const svgRef = useRef();
  const containerRef = useRef();
  const [dimensions, setDimensions] = useState({ width: initialWidth, height: height(initialWidth) });

  useEffect(() => {
    const updateDimensions = () => {
      if (containerRef.current) {
        const { width } = containerRef.current.getBoundingClientRect();
        setDimensions({ width, height: height(width) || 400 });
      }
    };

    updateDimensions();
    window.addEventListener("resize", updateDimensions);
    return () => window.removeEventListener("resize", updateDimensions);
  }, []);

  useEffect(() => {
    const svg = d3.select(svgRef.current);
    const { width, height } = dimensions;
    const margin = { top: 20, right: 30, bottom: 40, left: 50 };

    svg.selectAll("*").remove();

    const xScale = d3
      .scaleTime()
      .domain(d3.extent(timestamps))
      .range([margin.left, width - margin.right]);

    const yScale = d3
      .scaleLinear()
      .domain([d3.min(values.flat()), d3.max(values.flat())])
      .nice()
      .range([height - margin.bottom, margin.top]);

    const xAxis = d3.axisBottom(xScale).ticks(10);
    const yAxis = d3.axisLeft(yScale);

    svg
      .append("g")
      .attr("transform", `translate(0,${height - margin.bottom})`)
      .call(xAxis);
    svg
      .append("g")
      .attr("transform", `translate(${margin.left},0)`)
      .call(yAxis);

    values.forEach((v, i) => {
      const line = d3
        .line()
        .x((d, i) => xScale(timestamps[i]))
        .y((d) => yScale(d))
        .curve(d3.curveMonotoneX);

      svg
        .append("path")
        .attr("class", "graph")
        .datum(v)
        .attr("fill", "none")
        .attr("stroke", colors[i])
        .attr("stroke-width", 1.5)
        .attr("d", line);
    });

    // Append mouseLine group
    const mouseLineGroup = svg
      .append("g")
      .attr("class", "mouseLineGroup")
      .style("opacity", "0")
      .style("pointer-events", "none");

    mouseLineGroup
      .append("path")
      .attr("class", "mouseLine")
      .style("stroke", "white")
      .style("stroke-width", "1px");

    values.forEach((_, i) => {
      mouseLineGroup
        .append("circle")
        .attr("class", "mouseCircle")
        .attr("r", 5)
        .style("stroke", colors[i])
        .style("fill", "none")
        .style("stroke-width", "1px");
    });

    const tooltip = d3
      .select(svgRef.current.parentNode) // Append tooltip to the container, scoped to this graph
      .append("div")
      .style("position", "absolute")
      .style("background", "rgba(255,255,255,0.9)")
      .style("padding", "5px 10px")
      .style("border", "1px solid #ccc")
      .style("border-radius", "4px")
      .style("box-shadow", "0 0 5px rgba(0,0,0,0.2)")
      .style("display", "none")
      .style("pointer-events", "none")
      .style("font-size", "12px");

    const bisectDate = d3.bisector((d) => d).center;

    svg.on("mousemove", (event) => {
      const [mouseX] = d3.pointer(event);
      const hoveredX = xScale.invert(mouseX);

      const closestIndex = bisectDate(timestamps, hoveredX.getTime());
      const closestTimestamp = timestamps[closestIndex];
      const yValues = values.map((v) => v[closestIndex]);

      if (closestTimestamp !== undefined) {
        tooltip
          .style("display", "block")
          .style("color", "black")
          .style("left", `${event.pageX + 10}px`)
          .style("top", `${event.pageY - 20}px`)
          .html(
            `<strong>${new Date(closestTimestamp).toLocaleString()}</strong><br>
            <div style="text-align:left">${graphNames.map((name, i) => `<b style="color: ${colors[i]}">${name}</b>: ${yValues[i].toFixed(2)}`).join("<br>")}</div>
          `
          );

        mouseLineGroup.style("opacity", "1");

        const yRange = yScale.range();
        const xCoor = xScale(closestTimestamp);

        mouseLineGroup.select(".mouseLine").attr("d", `M${xCoor},${yRange[0]}L${xCoor},${yRange[1]}`);

        mouseLineGroup.selectAll(".mouseCircle").each(function (d, i) {
          const yCoor = yScale(yValues[i]);
          d3.select(this).attr("transform", `translate(${xCoor},${yCoor})`);
        });
      }
    });

    svg.on("mouseout", () => {
      mouseLineGroup.style("opacity", "0");
      tooltip.style("display", "none");
    });
  }, [timestamps, values, graphNames, dimensions]);

  return (<div ref={containerRef} style={{ width: "100%", height: "100%" }}>
    <svg
      ref={svgRef}
      width={dimensions.width}
      height={dimensions.height}
      viewBox={`0 0 ${dimensions.width} ${dimensions.height}`}
      preserveAspectRatio="xMinYMin meet"
    />
  </div>);
};

export default LineGraph;
