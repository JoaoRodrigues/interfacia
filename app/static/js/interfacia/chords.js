var datasets = [{
    "name": "van der Waals",
    "data": window.vdw_raw_data
}, {
    "name": "Electrostatics",
    "data": window.elec_raw_data
}]

// Setup SVG Canvas
var width = Math.round(parseInt(d3.select('#chart').style('width'), 10)),
    height = Math.round(parseInt(d3.select('#chart').style('height'), 10)),
    outerRadius = Math.min(width, height) / 2 - 30,
    innerRadius = outerRadius - 24;

// Create Plot
var layout = d3.layout.chord()
    .padding(.04)
    .sortSubgroups(d3.descending)
    .sortChords(d3.ascending);

var arc_svg = d3.svg.arc()
    .innerRadius(innerRadius)
    .outerRadius(outerRadius);

var path_svg = d3.svg.chord()
    .radius(innerRadius);

// Plotting function
function plot(nodes, raw_data) {

    // Process Input Data
    var o = scale_data(raw_data);
    var scaled_data = o[0];
    var mapping = o[1];
    var colors = get_color_scale(raw_data);

    // Get chain colors
    // TODO

    d3.select("#chart").selectAll('svg').remove(); // Clear 'canvas'

    var svg = d3.select("#chart").append("svg")
        .attr("width", width)
        .attr("height", height)
        .append("g")
        .attr("id", "circle")
        .attr("transform", "translate(" + width / 2 + "," + height / 2 + ")");

    svg.append("circle")
        .attr("r", outerRadius);

    // Compute the chord/group layout
    layout.matrix(scaled_data);

    // Add a group per interface residue
    var group = svg.selectAll(".group")
        .data(layout.groups)
        .enter().append("g")
        .attr("class", "group")
        .on("mouseover", mouseover)
        .on("click", group_click)


    var group_path = group.append("path")
        .attr("id", function(d, i) {
            return "group" + i;
        })
        .attr("d", arc_svg)
        // Fill
        .style("fill", function(d, i) {
            if (nodes[i].seg == 'A') {
                return "#006600";
            } else {
                return "#ff9900";
            }
        })
        // Mouseover Title
    group_path.append("title").text(function(d, i) {
        return nodes[i].seg + ":" + nodes[i].resi + ":" + nodes[i].resn;
    });

    // Add a chord per interaction
    var chord = svg.selectAll(".chord")
        .data(layout.chords)
        .enter().append("path")
        .attr("class", "chord")
        .attr("id", function(d, i) {
            return "chord" + i;
        })
        .style("fill", function(d) {
            return colors(mapping[d.source.value])
        })
        .attr("d", path_svg)
        .on("click", chord_click);

    // Add an elaborate mouseover title for each chord.
    chord.append("title").text(function(d) {
        return nodes[d.source.index].seg + ":" + nodes[d.source.index].resi + ":" + nodes[d.source.index].resn + " → " + nodes[d.target.index].seg + ":" + nodes[d.target.index].resi + ":" + nodes[d.target.index].resn + ": " + mapping[d.source.value];
    });

    // Mouseover function to fade chords
    function mouseover(d, i) {
        chord.classed("fade", function(p) {
            return p.source.index != i && p.target.index != i;
        });
    };

    // Clickers & Pass to PV selection and highlight
    function group_click(d, i) {
        // Single Residue
        var resViewName = 'pickedRes_' + nodes[i].seg + nodes[i].resi
        pickedRes = window.structure.select({
            chain: nodes[i].seg,
            rnum: nodes[i].resi
        })
        if (viewer.get(resViewName) === null) {
            viewer.ballsAndSticks(resViewName, pickedRes);
            d3.select("#group" + i)
                .classed("highlight", true);

        } else {
            // Hide Spheres
            viewer.rm(resViewName);
            // Highlight Group
            d3.select("#group" + i)
                .classed("highlight", false);
        }
        viewer.requestRedraw();
    };

    function chord_click(d, i) {
        // Bond (default to CA)
        var sourceAtom = window.structure.atom(nodes[d.source.index].seg + '.' + nodes[d.source.index].resi + '.CA')
        if (sourceAtom == null) {
            sourceAtom = window.structure.atom(nodes[d.source.index].seg + '.' + nodes[d.source.index].resi + '.C4')
        }

        var targetAtom = window.structure.atom(nodes[d.target.index].seg + '.' + nodes[d.target.index].resi + '.CA')
        if (targetAtom == null) {
            targetAtom = window.structure.atom(nodes[d.target.index].seg + '.' + nodes[d.target.index].resi + '.C4')
        }

        var bondViewName = 'pickedBond_' + sourceAtom.qualifiedName() + "&" + targetAtom.qualifiedName()
        if (viewer.get(bondViewName) === null) {
            var g = viewer.customMesh(bondViewName);
            var midPoint = pv.vec3.clone(sourceAtom.pos());
            pv.vec3.add(midPoint, midPoint, targetAtom.pos());
            pv.vec3.scale(midPoint, midPoint, 0.5);
            // add a tube to connect the two atoms
            // Color based on energy
            g.addTube(sourceAtom.pos(), targetAtom.pos(), 1.0,
                // { cap : true, color : colors(mapping[d.source.value]) });
                {
                    cap: true,
                    color: '#ffcc00'
                });

            d3.select('#chord' + i)
                .classed('highlight', true)
        } else {
            // Hide
            viewer.rm(bondViewName);
            d3.select('#chord' + i)
                .classed('highlight', false)
        }
        viewer.requestRedraw();
    };

    // Ticks & Labels
    group.append("text")
        .each(function(d) {
            d.angle = (d.startAngle + d.endAngle) / 2;
        })
        .attr("class", "chord_label")
        .attr("dy", ".35em")
        .attr("transform", function(d) {
            return "rotate(" + (d.angle * 180 / Math.PI - 90) + ")" + "translate(" + (innerRadius + 26) + ")" + (d.angle > Math.PI ? "rotate(180)" : "");
        })
        .style("text-anchor", function(d) {
            return d.angle > Math.PI ? "end" : null;
        })
        .text(function(d, i) {
            return nodes[i].resi;
        });

};

// Default to vdW data
plot(nodes, vdw_raw_data);

// Dataset switcher
var dataset_picker = d3.select("#dataset-picker").selectAll(".dataset-button")
    .data(datasets);

dataset_picker.enter()
    .append("input")
    .attr("value", function(d) {
        return "Dataset " + d.name
    })
    .attr("type", "button")
    .attr("class", "dataset-button")
    .on("click", function(d) {
        plot(nodes, d.data);
    });

// Transitions
// TODO

// Ramp
// TODO

// Function to scale raw data
// Scales each data point to the % of *total* *absolute* energy
// This means both highly favourable and unfavourable values will have similar percentages
// Sort out between them using color.
function scale_data(mtx) {
    var _abs_array = mtx.map(function(subArray) {
        return subArray.map(function(el) {
            return Math.abs(el)
        })
    })
    var _abs_total_ene = d3.sum(_abs_array.map(function(subArray) {
        return d3.sum(subArray)
    }))
    var _scaled_array = mtx.map(function(subArray) {
        return subArray.map(function(el) {
            return (Math.abs(el) / _abs_total_ene)
        })
    })

    var _flattened_raw_array = [].concat.apply([], mtx)
    var _flattened_scaled_array = [].concat.apply([], _scaled_array)
    var _mapping = []
    for (var i = 0; i < _flattened_raw_array.length; i++) {
        _mapping[_flattened_scaled_array[i]] = _flattened_raw_array[i]
    }

    return [_scaled_array, _mapping]
};

// Build color scale from energy values
function get_color_scale(mtx) {
    var _flattened_array = [].concat.apply([], mtx)
    var extremes = d3.extent(_flattened_array)
    var c = d3.scale.linear()
        .domain([extremes[0], 0, extremes[1]])
        // .range(["blue", "white", "red"]) // Flat
        .range(["#138BFF", "#fff", "#FF0007"]) // Vivid
        // .range(["#356089", "#fff", "#D1464A"]) // Pastel
    return c
};
