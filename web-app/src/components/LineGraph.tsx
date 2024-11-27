import React, { useEffect, useRef } from "react";
import * as d3 from "d3";

const colors = ["steelblue", "orange", "green", "#ff6361"];

const LineGraph = ({ timestamps, values, graphNames }) => {
  const svgRef = useRef();

  useEffect(() => {
    const svg = d3.select(svgRef.current);
    const width = 800;
    const height = 400;
    const margin = { top: 20, right: 30, bottom: 40, left: 50 };

    svg.selectAll("*").remove();

    const xScale = d3
      .scaleTime()
      .domain(d3.extent(timestamps))
      .range([margin.left, width - margin.right]);

    const yScale = d3
      .scaleLinear()
      .domain([
        d3.min(values.flat()),
        d3.max(values.flat()),
      ])
      .nice()
      .range([height - margin.bottom, margin.top]);

    const xAxis = d3.axisBottom(xScale).ticks(10);
    const yAxis = d3.axisLeft(yScale);

    svg.append("g").attr("transform", `translate(0,${height - margin.bottom})`).call(xAxis);
    svg.append("g").attr("transform", `translate(${margin.left},0)`).call(yAxis);

    

    const line1 = d3
      .line()
      .x((d, i) => xScale(timestamps[i]))
      .y((d) => yScale(d))
      .curve(d3.curveMonotoneX);

    const line2 = d3
      .line()
      .x((d, i) => xScale(timestamps[i]))
      .y((d) => yScale(d))
      .curve(d3.curveMonotoneX);

    const line3 = d3
      .line()
      .x((d, i) => xScale(timestamps[i]))
      .y((d) => yScale(d))
      .curve(d3.curveMonotoneX);

    values.forEach((v, i) => {
      svg.append("path").attr('class', 'graph').datum(v).attr("fill", "none").attr("stroke", colors[i]).attr("stroke-width", 1.5).attr("d", line1);

    });

    const mouseLine = svg.append('g').attr("class","mouseLineGroup").style("opacity", "0");

    mouseLine
      .append("path")
      .attr("class","mouseLine")
      .style("stroke","white")
      .style("stroke-width", "1px");

    svg.selectAll('.graph').each((_, i) => {
      mouseLine.append("circle")
      .attr("class","mouseCircle") // add a circle to follow along path
      .attr("r", 5)
      .style("stroke", colors[i])
      .style("fill","none")
      .style("stroke-width", "1px")
    })

    const tooltip = d3
      .select("body")
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

    const overlay = svg
      .append("rect")
      .attr("fill", "none")
      .attr("pointer-events", "all")
      .attr("width", width - margin.left - margin.right)
      .attr("height", height - margin.top - margin.bottom)
      .attr("transform", `translate(${margin.left},${margin.top})`);

    const bisectDate = d3.bisector((d) => d).center;

    overlay.on("mousemove", (event) => {
      const [mouseX] = d3.pointer(event);
      const hoveredX = xScale.invert(mouseX + margin.left);
      
      const closestIndex = bisectDate(timestamps, hoveredX.getTime());
      
      const closestTimestamp = timestamps[closestIndex];
      const yValues = values.map(v => v[closestIndex]);

      if (closestTimestamp !== undefined) {
        tooltip
          .style("display", "block")
          .style("color", "black")
          .style("left", `${xScale(closestTimestamp) + 10}px`)
          .style("top", `${event.pageY - 20}px`)
          .html(`
            <strong>${new Date(closestTimestamp).toLocaleString()}</strong><br>
            ${graphNames[0]}: ${yValues[0].toFixed(2)}<br>
            ${graphNames[1]}: ${yValues[1].toFixed(2)}<br>
            ${graphNames[2]}: ${yValues[2].toFixed(2)}
          `);
      }

      d3.select(".mouseLine")
      .attr("d", function(){
          const yRange = yScale.range(); // range of y axis
          const xCoor = xScale(closestTimestamp); // mouse position in x
          
          d3.selectAll('.mouseCircle') // for each circle group
              .each(function(d,i){
                const yCoor = yScale(yValues[i]);
                d3.select(this) // move the circle to intersection
                  .attr('transform', 'translate(' + xCoor + ',' + yCoor + ')');
              });
          return "M"+ xCoor +"," + yRange[0] + "L" + xCoor + "," + yRange[1]; // position vertical line
      });
    });

    overlay.on('mouseover', function(){ // on mouse in show line, circles and text
      d3.select(".mouseLineGroup")
          .style("opacity", "1");
    }).on("mouseout", () => {
      d3.select(".mouseLineGroup")
        .style("opacity", "0");
      tooltip.style("display", "none");
    });
  }, [timestamps, values, graphNames]);

  return <svg ref={svgRef} width={800} height={400}></svg>;
};

export default LineGraph;
