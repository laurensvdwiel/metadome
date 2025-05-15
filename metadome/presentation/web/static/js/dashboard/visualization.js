/*******************************************************************************
 * Global variables for landscape
 ******************************************************************************/
var selected_positions = 0;
var meta_domain_ids = new Set();

var main_outerWidth = 1300;
var main_outerHeight = 500;
var main_svg = d3.select("#landscape_svg").attr("width", main_outerWidth)
		.attr("height", main_outerHeight);

// Declare margins
var main_marginLandscape = {
	top : 20,
	right : 20,
	bottom : 210,
	left : 80
};
var main_marginLegend = {
	top : 20,
	right : 1240,
	bottom : 210,
	left : 20
};
var main_marginPositionInfo = {
	top : 325,
	right : 20,
	bottom : 135,
	left : 80
};
var main_marginAnnotations = {
	top : 310,
	right : 20,
	bottom : 120,
	left : 80
};
var main_marginContext = {
	top : 410,
	right : 20,
	bottom : 30,
	left : 80
};

// Declare various UI widths and heights
var main_width = main_outerWidth - main_marginLandscape.left - main_marginLandscape.right;
var main_widthLegend = main_outerWidth - main_marginLegend.left - main_marginLegend.right;
var main_heightLandscape = main_outerHeight - main_marginLandscape.top
		- main_marginLandscape.bottom;
var main_heightMarginPositionInfo = main_outerHeight - main_marginPositionInfo.top
		- main_marginPositionInfo.bottom;
var main_heightContext = main_outerHeight - main_marginContext.top - main_marginContext.bottom;
var main_heightAnnotations = main_outerHeight - main_marginAnnotations.top
		- main_marginAnnotations.bottom;

// Scale the axis
var main_x = d3.scaleLinear().range([ 40, main_width -40 ]).nice();
var main_x2 = d3.scaleLinear().range([ 0, main_width ]);
var main_y = d3.scaleLinear().range([ main_heightLandscape, 0 ]);
var main_y_metadomain = d3.scaleLinear().range([ main_heightLandscape, 0 ]);

// Set axis
var main_xAxis = d3.axisBottom(main_x2).ticks(0);
var main_yAxis = d3.axisLeft(main_y).ticks(0);

// add the brush element
var brush = d3.brushX().extent([ [ 0, 0 ], [ main_width, main_heightContext ] ]).on("brush end", brushed);
// Define the brushed function
function brushed() {
	var s = d3.event.selection || main_x2.range();
	main_x.domain(s.map(main_x2.invert, main_x2));
	rescaleLandscape();
}
/*******************************************************************************
 * Config variables for visuals
 ******************************************************************************/

// indicates the maximum tolerance score
var maxTolerance = 1.8;

// indicates if the metadomain landscape is visible
var metadomain_graph_visible = true;

// indicates if clinvar variants are annotated in the schematic protein representation
var clinvar_variants_visible = false;

//indicates if homologous clinvar variants are annotated in the schematic protein representation
var homologous_clinvar_variants_visible = false;

// indicates the various colors to indicate the tolerance
var toleranceColorGradient = [ {
	offset : "0%",
	color : "#d7191c"
}, {
	offset : "12.5%",
	color : "#e76818"
}, {
	offset : "25%",
	color : "#f29e2e"
}, {
	offset : "37.5%",
	color : "#f9d057"
}, {
	offset : "50%",
	color : "#ffff8c"
}, {
	offset : "62.5%",
	color : "#90eb9d"
}, {
	offset : "75%",
	color : "#00ccbc"
}, {
	offset : "87.5%",
	color : "#00a6ca"
}, {
	offset : "100%",
	color : "#2c7bb6"
} ]

/*******************************************************************************
 * Basic user interface elements
 ******************************************************************************/

// Define Area of Tolerance landscape graph
var toleranceArea = d3.area().x(function(d) {
	return main_x (d.protein_pos);
}).y0(main_heightLandscape).y1(function(d) {
	return main_y(d.sw_dn_ds);
});

// Define Line of Tolerance landscape line plot
var toleranceLine = d3.line().x(function(d) {
	return main_x (d.values[0].protein_pos);
}).y(function(d) {
	return main_y(d.values[0].sw_dn_ds);
});

/*******************************************************************************
 * Additional user interface elements
 ******************************************************************************/

// Define tooltip for pfamdomains
var domainTip = d3.tip()
	.attr('class', 'd3-tip')
	.offset([ -10, 0 ])
	.html(function(d) {
	    return "<span><div style='text-align: center;'> Domain " + d.Name + " ("+d.ID+") </br> click for more information </div></span>";
	});

//Define tooltip for pfamdomains
var domain_details_position_tip = d3.tip()
	.attr('class', 'd3-tip')
	.offset([ -10, 0 ])
	.html(function(d) {
	    return "<span> "+ d +"</span>";
	});


// Define tooltip for positions
var positionTip = d3.tip()
	.attr('class', 'd3-tip')
	.offset([ -10, 0 ])
	.html(function(d, i) {
	    var positionTip_str = "<span>";
	    positionTip_str += "Position: p." + d.values[0].protein_pos + " " + d.values[0].cdna_pos + "</br>";
	    positionTip_str += "Codon: " + d.values[0].ref_codon + "</br>";
	    positionTip_str += "Residue: " + d.values[0].ref_aa_triplet + "</br>";
	    positionTip_str += "Tolerance score (dn/ds): "+ (Math.round((d.values[0].sw_dn_ds)*100)/100) +' ('+tolerance_rating(d.values[0].sw_dn_ds) +')';
	    if (d.values[0].domains.length > 0){
		positionTip_str += "</br> In domain(s): ";
		var n_domains_at_position = d.values[0].domains.length;
		for (var i = 0; i < n_domains_at_position; i++){
		    if (i+1 == n_domains_at_position){
			positionTip_str+= d.values[0].domains[i].ID;
		    }else{
			positionTip_str+= d.values[0].domains[i].ID+", ";
		    }
		}
	    }
	    positionTip_str += "</br> Click to select this position</span>";
	    return positionTip_str;
	});

/*******************************************************************************
 * Main draw function
 ******************************************************************************/

//Reset all graph elements based on the obj
function resetGraph(){
    // reset the variables
	selected_positions = 0;
	meta_domain_ids = new Set();
	
	$("#selected_positions_information").addClass('is-hidden');
	document.getElementById("selected_positions_explanation").innerHTML = 'Click on one of positions in the schematic protein to obtain more information';
	$("#position_information_table").addClass('is-hidden');
	d3.selectAll('.tr').remove();

	// reset the svg
	main_svg.selectAll("*").remove();
	main_svg = d3.select("#landscape_svg").attr("width", main_outerWidth)
	.attr("height", main_outerHeight);
}

// Creates all graph elements based on the obj
function createGraph(obj) {
	$("#geneName").html(obj.geneName);
	
	// reset the Graph
	resetGraph();
	
	// Make the Selected positions visible
	$("#selected_positions_information").removeClass('is-hidden');
	
	// Add defs to the svg
	var defs = main_svg.append("defs")
		.attr('class', 'defs');

	// Add clipping to defs
	main_svg.select('.defs').append("clipPath")
		.attr("id", "clip")
		.append("rect")
		.attr("width", main_width)
		.attr("height", main_heightLandscape);

	// Extract the various data
	var positional_annotation = obj.positional_annotation;
	var domain_data = obj.domains;

	// setting x/y domain according to data
	main_x.domain(d3.extent(positional_annotation, function(d) {
	    return d.protein_pos;
	}));
	main_y.domain([ 0, maxTolerance ]);
	main_x2.domain(main_x .domain());

	// create a group from the tolerance data
	var dataGroup = d3.nest().key(function(d) {
	    return d.protein_pos;
	}).entries(positional_annotation);

	// add two consecutive data values per group, so these can be used in
	// drawing rectangles
	dataGroup.forEach(function(group, i) {
	    	group.values.push([]);
		if (i + 1 < dataGroup.length) {
		    group.values[1].protein_pos = dataGroup[i + 1].values[0].protein_pos;
		    group.values[1].sw_dn_ds = dataGroup[i + 1].values[0].sw_dn_ds;
		} else if (i < dataGroup.length) {
		    group.values[1].protein_pos = dataGroup[i].values[0].protein_pos;
		    group.values[1].sw_dn_ds = dataGroup[i].values[0].sw_dn_ds;
		}
		group.values[0].selected = false;
	})
	
	// create an overview of coverage per domain ID
	var domain_metadomain_coverage = {}
	for (var i = 0; i < domain_data.length; i++) { 
		if (!(domain_data[i].ID in domain_metadomain_coverage)){
			domain_metadomain_coverage[domain_data[i].ID] = domain_data[i].meta_domain_alignment_depth;
		}
	}
	
	// Draw all individual user interface elements based on the data
	annotateDomains(domain_data, positional_annotation, domain_metadomain_coverage);
	createToleranceGraph(dataGroup);
	createToleranceGraphLegend();
	drawMetaDomainLandscape(domain_data, dataGroup, domain_metadomain_coverage, obj.transcript_id);
	createMetaDomainLegend();

	// Add schematic protein overview as a custom Axis
	createSchematicProtein(domain_metadomain_coverage, dataGroup, obj.transcript_id);

	// Finally draw the context zoom
	addContextZoomView(domain_data, dataGroup.length);
	
	// Add behaviour according to the settings
	toggleToleranceLandscapeOrMetadomainLandscape();
}

function drawMetaDomainLandscape(domain_data, data, domain_metadomain_coverage, transcript_id){
    // get all possible domain ids
    for (var i = 0; i < domain_data.length; i++){
    	if (domain_data[i].metadomain){
    		meta_domain_ids.add(domain_data[i].ID);
    	}
    }
    
    // receive the max pathogenic and normal variation
    var global_max_normal = 0;
    var global_max_pathogenic = 0;
    var global_max_value = 0;
    for (var i = 0; i < data.length; i++){
    	meta_domain_ids.forEach(domain_id => {
			if (data[i].values[0].hasOwnProperty('domains') && data[i].values[0].domains[domain_id] != null){
			    global_max_normal = Math.max(data[i].values[0].domains[domain_id].normal_missense_variant_count, global_max_normal);
			    global_max_pathogenic = Math.max(data[i].values[0].domains[domain_id].pathogenic_missense_variant_count, global_max_pathogenic);
			}
		});    	
    }
	if (global_max_pathogenic == 0){
		global_max_value = global_max_normal;
	}
	else {
		global_max_value = global_max_pathogenic;
	}
        
	// Add barplot for the metadomain variation landscape
	var meta_domain_landscape_canvas = main_svg.append("g")
    .attr("id", "metadomain_graph")
    .attr("transform", "translate(" + main_marginLandscape.left + "," + main_marginLandscape.top + ")");
	
	// Call the tooltips
    meta_domain_landscape_canvas.call(domain_details_position_tip);

    
	// Define the axes based on the data
	main_y_metadomain.domain([0, global_max_value]);
	
	// Draw the meta-domain normal missense variation barplot
	meta_domain_landscape_canvas.selectAll(".bar")
	.data(data)
	.enter()
	.append("rect")
	.attr("class", "normal_missense_variant_count")
	.attr("x", function(d) { return main_x(d.values[0].protein_pos - 0.5); })
	.attr("y", function(d) { 
		var normal_missense_variant_count = 0;
		if (d.values[0].domains != null){
			meta_domain_ids.forEach(domain_id => {
				if (d.values[0].hasOwnProperty('domains') && d.values[0].domains[domain_id] != null){
					normal_missense_variant_count = Math.max(d.values[0].domains[domain_id].normal_missense_variant_count, normal_missense_variant_count);
				}
			});

			return main_y_metadomain((global_max_value/global_max_normal)*normal_missense_variant_count);
		}
		})
	.attr("width", function(d, i) {
	    if (d.values[1].protein_pos != d.values[0].protein_pos){
			return main_x(d.values[1].protein_pos) - main_x(d.values[0].protein_pos);
	    	    } else {
	    		return main_x(d.values[1].protein_pos +1) - main_x(d.values[0].protein_pos);
	    	    }
		})
	.attr("height", function(d) { 
		var normal_missense_variant_count = 0;
		if (d.values[0].domains != null){
			meta_domain_ids.forEach(domain_id => {
				if (d.values[0].hasOwnProperty('domains') && d.values[0].domains[domain_id] != null){
					normal_missense_variant_count = Math.max(d.values[0].domains[domain_id].normal_missense_variant_count, normal_missense_variant_count);
				}
			});

			return main_heightLandscape - main_y_metadomain((global_max_value/global_max_normal)*normal_missense_variant_count);
		}
		})
	.style("clip-path", "url(#clip)")
	.style("fill", "green")
	.on("click", function(d) {
	    // Call this method found in dashboard.js
	    createPositionalInformation(domain_metadomain_coverage, transcript_id, d)
	})
	.on("mouseover", function(d) {
		if (metadomain_graph_visible){
			var normal_missense_variant_count = 0;
			if (d.values[0].domains != null){
				meta_domain_ids.forEach(domain_id => {
					if (d.values[0].hasOwnProperty('domains') && d.values[0].domains[domain_id] != null){
						normal_missense_variant_count = Math.max(d.values[0].domains[domain_id].normal_missense_variant_count, normal_missense_variant_count);
					}
				});
			}
		   // show the tooltip
		   domain_details_position_tip.show("Homologous gnomAD missense count: "+normal_missense_variant_count);
		   // amplify the element
		   d3.select(this).style("fill", "orange");
		   // move the element to front
		   d3.select(this).moveToFront();
		}
	})
	.on("mouseout", function(d) {
	   // hide the tooltip
	   domain_details_position_tip.hide(d);
	   // reset the color
	   d3.select(this).style("fill", "green");
	   // move the element to the back
	   d3.select(this).moveToBack();
	});
	
	// Draw the meta-domain pathogenic missense variation barplot
	meta_domain_landscape_canvas.selectAll(".bar")
	.data(data)
	.enter()
	.append("rect")
	.attr("class", "pathogenic_missense_variant_count")
	.attr("x", function(d) { return main_x(d.values[0].protein_pos); })
	.attr("y", function(d) { 
		var pathogenic_missense_variant_count = 0;
		if (d.values[0].domains != null){
			meta_domain_ids.forEach(domain_id => {
				if (d.values[0].hasOwnProperty('domains') && d.values[0].domains[domain_id] != null){
					pathogenic_missense_variant_count = Math.max(d.values[0].domains[domain_id].pathogenic_missense_variant_count, pathogenic_missense_variant_count);
				}
			});

			return main_y_metadomain((global_max_value/global_max_pathogenic)*pathogenic_missense_variant_count);
		}
		})
	.attr("width", function(d, i) {
	    if (d.values[1].protein_pos != d.values[0].protein_pos){
			return main_x(d.values[1].protein_pos) - main_x(d.values[0].protein_pos);
	    	    } else {
	    		return main_x(d.values[1].protein_pos +1) - main_x(d.values[0].protein_pos);
	    	    }
		})
	.attr("height", function(d) { 
		var pathogenic_missense_variant_count = 0;
		if (d.values[0].domains != null){
			meta_domain_ids.forEach(domain_id => {
				if (d.values[0].hasOwnProperty('domains') && d.values[0].domains[domain_id] != null){
					pathogenic_missense_variant_count = Math.max(d.values[0].domains[domain_id].pathogenic_missense_variant_count, pathogenic_missense_variant_count);
				}
			});
			
			return main_heightLandscape - main_y_metadomain((global_max_value/global_max_pathogenic)*pathogenic_missense_variant_count);
		}
		})
	.style("clip-path", "url(#clip)")
	.style("fill", "red")
	.on("click", function(d) {
	    // Call this method found in dashboard.js
	    createPositionalInformation(domain_metadomain_coverage, transcript_id, d)
	})
	.on("mouseover", function(d) {
		if (metadomain_graph_visible){
			var pathogenic_missense_variant_count = 0;
			if (d.values[0].domains != null){
				meta_domain_ids.forEach(domain_id => {
					if (d.values[0].hasOwnProperty('domains') && d.values[0].domains[domain_id] != null){
						pathogenic_missense_variant_count = Math.max(d.values[0].domains[domain_id].pathogenic_missense_variant_count, pathogenic_missense_variant_count);
					}
				});
			}
		   // show the tooltip
		   domain_details_position_tip.show("Homologous pathogenic missense count: "+pathogenic_missense_variant_count);
		   // amplify the element
		   d3.select(this).style("fill", "orange");
		   // move the element to front
		   d3.select(this).moveToFront();
		}
	})
	.on("mouseout", function(d) {
	   // hide the tooltip
	   domain_details_position_tip.hide(d);
	   // reset the color
	   d3.select(this).style("fill", "red");
	   // move the element to the back
	   d3.select(this).moveToBack();
	});
	
	// Draw not aligned barplots
	meta_domain_landscape_canvas.selectAll(".bar")
	.data(data)
	.enter()
	.append("rect")
	.attr("class", "not_aligned_position_plot")
	.attr("x", function(d) { return main_x(d.values[0].protein_pos); })
	.attr("y", function(d) { 
		if (d.values[0].domains != null){
			var not_aligned_poition = 0;
			
			for (var i = 0; i < domain_data.length; i++){
		    	if (domain_data[i].start <= d.values[0].protein_pos && d.values[0].protein_pos < domain_data[i].stop && d.values[0].domains[domain_data[i].ID] == null){
		    		not_aligned_poition = 1;
		    	}
		    }
			
			return main_y_metadomain((global_max_value/10)*not_aligned_poition);
		}
		})
	.attr("width", function(d, i) {
	    if (d.values[1].protein_pos != d.values[0].protein_pos){
			return main_x(d.values[1].protein_pos) - main_x(d.values[0].protein_pos);
	    	    } else {
	    		return main_x(d.values[1].protein_pos +1) - main_x(d.values[0].protein_pos);
	    	    }
		})
	.attr("height", function(d) { 
		if (d.values[0].domains != null){
			
			return main_heightLandscape - main_y_metadomain(global_max_value/10);
		}
		})
	.style("clip-path", "url(#clip)")
	.style("fill", "black")
	.on("mouseover", function(d) {
		if (metadomain_graph_visible){
			var pathogenic_missense_variant_count = 0;
			if (d.values[0].domains != null){
				meta_domain_ids.forEach(domain_id => {
					if (d.values[0].hasOwnProperty('domains') && d.values[0].domains[domain_id] != null){
						pathogenic_missense_variant_count = Math.max(d.values[0].domains[domain_id].pathogenic_missense_variant_count, pathogenic_missense_variant_count);
					}
				});
			}
		   // show the tooltip
		   domain_details_position_tip.show("Residue is not aligned to homologues");
		   // amplify the element
		   d3.select(this).style("fill", "orange");
		   // move the element to front
		   d3.select(this).moveToFront();
		}
	})
	.on("mouseout", function(d) {
	   // hide the tooltip
	   domain_details_position_tip.hide(d);
	   // reset the color
	   d3.select(this).style("fill", "black");
	   // move the element to the back
	   d3.select(this).moveToBack();
	});
}

/*******************************************************************************
 * Drawing user interface component
 ******************************************************************************/

// Draw the tolerance graph
function createToleranceGraph(dataGroup) {
	// append focus view
	var focus = main_svg.append("g")
		.attr("class", "focus")
		.attr("id", "tolerance_graph")
		.attr("transform", "translate(" + main_marginLandscape.left + "," + main_marginLandscape.top + ")");
	
	// draw the tolerance area graph, base on the grouped consecutive data
	// values
	dataGroup.forEach(function(d) {
	    // add the area specific for this position
	    focus.append("path").datum(d.values)
	    	.attr("class", "area")
	    	.attr("d", toleranceArea);

	    // create a linear gradient, specific for this position
	    var lineargradient = main_svg.append("linearGradient")
	    	.attr("id", "area-gradient_" + d.values[0].protein_pos + "-" + d.values[1].protein_pos)
	    	.attr("x1", "0%")
		.attr("y1", "0%")
		.attr("x2", "100%")
		.attr("y2", "0%");

	    // calculate the offset start score
	    lineargradient.append("stop")
	    	.attr("class", "start")
	    	.attr("offset", "0%")
	    	.attr("stop-color", tolerance_color(d.values[0].sw_dn_ds));

	    // calculate the offset stop score
	    lineargradient.append("stop")
	    	.attr("class", "end")
	    	.attr("offset", "100%")
	    	.attr("stop-color", tolerance_color(d.values[1].sw_dn_ds));
	});

	// color the area under the curve in contrast to the tolerance score
	focus.selectAll(".area")
		.style("fill", function(d, i) {
		    return "url(#area-gradient_" + d[0].protein_pos + "-" + d[1].protein_pos + ")";
		});

	// add tolerance line
	focus.append('path').datum(dataGroup)
		.attr('class', 'line')
		.attr('id', 'toleranceLine')
		.attr('fill', 'none')
		.attr("stroke", "black")
		.attr('stroke-width', "0.5px")
		.style("clip-path", "url(#clip)")
		.attr('d', toleranceLine);

	// append yAxis for focus view
	focus.append("g").attr("class", "axis axis--y").call(main_yAxis);
}

// Draw the axis and labels
function createSchematicProtein(domain_metadomain_coverage, groupedTolerance, transcript_id) {
	// Add the Axis
	var focusAxis = main_svg.append("g")
		.attr("class", "focusAxis")
		.attr("id", "tolerance_axis");

	// Add the elements that will be drawn
	var focusAxiselements = focusAxis.selectAll("g")
		.data(groupedTolerance).enter()
		.append("g")
		.attr("class", "focusAxisElement")
		.attr("transform", "translate(" + main_marginPositionInfo.left + ","	+ main_marginPositionInfo.top + ")")
		.style("fill", "none");
	
	// Call the tooltips
	main_svg.call(positionTip);

	focusAxiselements.append("rect")
	.attr("class", "toleranceAxisTickBackground")
	.attr("x", function(d, i) {
	    return main_x(d.values[0].protein_pos - 0.5);
	}).attr("y", 0).attr("width", function(d, i) {
	    if (d.values[1].protein_pos != d.values[0].protein_pos){
		return main_x(d.values[1].protein_pos) - main_x(d.values[0].protein_pos);
	    } else {
		return main_x(d.values[1].protein_pos +1) - main_x(d.values[0].protein_pos);
	    }
	})
	.attr("height", main_heightMarginPositionInfo)
	.style("fill", "white")
	.style("clip-path","url(#clip)");
	
	// Add a text per position
	focusAxiselements.append("text")
		.attr("class", "toleranceAxisTickLabel")
		.attr("id", function(d, i) {
		    return "toleranceAxisText_" + d.values[0].protein_pos;
		})
		.attr("x", function(d, i) {
		    return main_x(d.values[0].protein_pos);
		})
		.attr("y", main_heightMarginPositionInfo / 2)
		.attr("dy", ".35em")
		.style('pointer-events', 'none')
		.style('user-select', 'none')
		.attr("text-anchor", "middle")
		.style("fill", "black")
		.style("clip-path", "url(#clip)")
		.text(function(d, i) {
			return d.values[0].ref_aa;
		});

	// Add a rectangle per position
	focusAxiselements.append("rect")
		.attr("class", "toleranceAxisTick")
		.attr("id", function(d, i) {
		    return "toleranceAxisRect_" + d.values[0].protein_pos;
		})
		.attr("x", function(d, i) {
		    return main_x(d.values[0].protein_pos - 0.5);
		}).attr("y", 0).attr("width", function(d, i) {
		    if (d.values[1].protein_pos != d.values[0].protein_pos){
			return main_x(d.values[1].protein_pos) - main_x(d.values[0].protein_pos);
		    } else {
			return main_x(d.values[1].protein_pos +1) - main_x(d.values[0].protein_pos);
		    }
		})
		.attr("height", main_heightMarginPositionInfo)
		.style("fill-opacity", 0.2)
		.style("fill", "grey")
		.style("clip-path","url(#clip)")
		.on("mouseover", function(d, i) {
		    if (!d.values[0].selected) {
		    	d3.select(this).style("fill", "orange").style("fill-opacity", 0.5);
		    }
		    // show the tooltip
		    positionTip.show(d);
		    // change the cursor
		    d3.select(this).style("cursor", "pointer");
		}).on("mouseout", function(d, i) {			
		    if (!d.values[0].selected) {
				d3.select(this).style("fill", draw_position_schematic_protein(d, this));
		    }
		    // hide the tooltip
		    positionTip.hide(d);
		}).on("click", function(d, i) {
		    if (!d.values[0].selected) {
				d3.select(this).style("fill", "green").style("fill-opacity", 0.7);
				d.values[0].selected = true;
				
				addRowToPositionalInformationTable(domain_metadomain_coverage, d, transcript_id);
				selected_positions += 1;
				$("#position_information_table").removeClass('is-hidden');
			    document.getElementById("selected_positions_explanation").innerHTML = 'Click on one of the selected positions in the table to view more information';
		    } else {
				d3.select(this).style("fill", "orange").style("fill-opacity", 0.5);
				d.values[0].selected = false;
				d3.select("#positional_table_info_" + d.values[0].protein_pos).remove();
	
				selected_positions -= 1;
				if (selected_positions <= 0){
				    $("#position_information_table").addClass('is-hidden');
					document.getElementById("selected_positions_explanation").innerHTML = 'Click on one of positions in the schematic protein to obtain more information';
				}
		    }
		});
}

// Draw the domain annotation
function annotateDomains(protDomain, tolerance_data, domain_metadomain_coverage) {
	// append domain view
	var domains = main_svg.append("g")
		.attr("class", "domains")
		.attr("id", "domain_annotation")
		.attr("transform", "translate(" + main_marginAnnotations.left + "," + main_marginAnnotations.top + ")");

	// Add text to the ui element
	main_svg.append("text")
		.attr("text-anchor", "left")
		.attr("x", 0).attr("y", main_marginAnnotations.top + (main_heightAnnotations*(3/5))).attr("dy", 0)
		.attr("class", "label")
		.text("Protein")
		.style('pointer-events', 'none')
		.style('user-select', 'none');

	// Call the tooltips
	main_svg.call(domainTip);

	// Fill the ui element
	domains.selectAll(".rect")
		.data(protDomain).enter()
		.append("rect").attr("class", "pfamDomains")
		.attr("x", function(d) {
		    return main_x(d.start - 0.5);
		})
		.attr("y", 0)
		.attr("width", function(d) {
		    return main_x(d.stop + 1) - main_x(d.start);
		})
		.attr("height", main_heightAnnotations)
		.attr("rx", 10)
		.attr("ry", 10)
		.style('opacity', 0.5)
		.style('fill', '#c014e2')
		.style('stroke', 'black')
		.style("clip-path", "url(#clip)")
		.on("mouseover", function(d) {
		    // show the tooltip
		    domainTip.show(d);
		    // amplify the element
		    d3.select(this).style("fill", "yellow");
		    // move the element to front
		    d3.select(this).moveToFront();
		    // change the cursor
		    d3.select(this).style("cursor", "pointer");
		})
		.on("mouseout", function(d) {
		    // hide the tooltip
		    domainTip.hide(d);
		    // reset the color
		    d3.select(this).style("fill", "#c014e2");
		    // move the element to the back
		    d3.select(this).moveToBack();
		    // change the cursor
		    d3.select(this).style("cursor", "default");
		})
		.on("click", function(d) {
			// Activate the overlay
		    $("#domain_information_overlay").addClass('is-active');
		    
		    domain_information_overlay_content
		    
		    // Add title to box
		    document.getElementById("domain_information_overlay_title").innerHTML = '';
		    document.getElementById("domain_information_overlay_title").innerHTML += '<label class="label" >';
		    document.getElementById("domain_information_overlay_title").innerHTML += d.Name+' (' + d.ID + '), p.'+d.start+' - p.'+d.stop;
		    document.getElementById("domain_information_overlay_title").innerHTML += '</label>';

		    // Format the HTML in the correct format
		    document.getElementById("domain_information_overlay_content").innerHTML = '';
		    // Add domain occurrence details
		    document.getElementById("domain_information_overlay_content").innerHTML += '<label class="label" >';
		    document.getElementById("domain_information_overlay_content").innerHTML += 'Domain: ' + d.Name+' (<a href="https://www.ebi.ac.uk/interpro/entry/pfam/' + d.ID + '" target="_blank">' + d.ID + '</a>), located at p.'+d.start+' - p.'+d.stop;
		    document.getElementById("domain_information_overlay_content").innerHTML += ' in ';
		    document.getElementById("domain_information_overlay_content").innerHTML += document.getElementById("geneDetails").innerHTML +'.<br>';
		    
		    // if there is any, add meta-domain details
		    if (d.metadomain){
			    document.getElementById("domain_information_overlay_content").innerHTML += 'This domain has '+domain_metadomain_coverage[d.ID]+' homologous occurrences throughout the human genome.';
		    }
		    document.getElementById("domain_information_overlay_content").innerHTML += '</label>';

		});

	// function to move item to front of svg
	d3.selection.prototype.moveToFront = function() {
		return this.each(function() {
			this.parentNode.appendChild(this);
		});
	};

	// function to move item to back of svg
	d3.selection.prototype.moveToBack = function() {
		return this.each(function() {
			var firstChild = this.parentNode.firstChild;
			if (firstChild) {
				this.parentNode.insertBefore(this, firstChild);
			}
		});
	};
}

// Draw the context area for zooming
function addContextZoomView(domain_data, number_of_positions) {
	
	// append context view
	var context = main_svg.append("g")
		.attr("class", "context")
		.attr("id", "zoom_landscape")
		.attr("transform", "translate(" + main_marginContext.left + "," + main_marginContext.top + ")");

	// append context area
	context.append("rect")
		.attr("x", function(d) {
		    return main_x2(0);
		})
		.attr("y", 20)
		.attr("width", function(d) {
		    return main_x2(number_of_positions);
		})
		.attr("height", main_heightContext-40)
		.attr("rx", 10)
		.attr("ry", 10)
		.style("fill", "grey")
		.style("clip-path", "url(#clip)");
	
	// Append each domain
	context.selectAll(".rect")
		.data(domain_data).enter()
		.append("rect")
		.attr("x", function(d) {
		    return main_x2(d.start - 0.5);
		})
		.attr("y", 10)
		.attr("width", function(d) {
		    return main_x2(d.stop + 1) - main_x2(d.start);
		})
		.attr("height", main_heightContext-20)
		.attr("rx", 10)
		.attr("ry", 10)
		.style("fill", "grey")
		.style("clip-path", "url(#clip)");

	// append xAxis for context view
	context.append("g")
		.attr("class", "axis axis--x")
		.attr("transform", "translate(0," + main_heightContext + main_marginContext.top + ")")
		.call(main_xAxis);

	// append yAxis for context view
	context.append("g")
		.attr("class", "brush")
		.attr("id", "brush_for_zooming")
		.call(brush)
		.call(brush.move);

	main_svg.append("text")
		.attr("text-anchor", "left")
		.attr("id", "schematic_protein_zoom_text")
		.attr("x", 0)
		.attr("y", main_marginContext.top + (main_heightContext*(3/5)))
		.attr("dy", 0)
		.attr("class", "label")
		.text("Zoom-in")
		.style('pointer-events', 'none')
		.style('user-select', 'none');

	
}

// Draw the legend for the Tolerance graph
function createToleranceGraphLegend() {
	// append gradient to defs
	var legendGradient = main_svg.select('.defs')
		.append("linearGradient")
		.attr("id", "legendGradient")
		.attr("x1", "0%")
		.attr("y1", "100%")
		.attr("x2", "0%")
		.attr("y2", "0%");

	// set legend of the tolerance score
	legendGradient.selectAll("stop").data(toleranceColorGradient).enter()
		.append("stop")
		.attr("offset", function(d) {
		    return d.offset;
		})
		.attr("stop-color", function(d) {
		    return d.color;
		});

	// append heatmap legend
	main_svg.append("rect")
		.attr("id", "legendGradientRect")
		.attr("width", main_widthLegend)
		.attr("height", main_heightLandscape)
		.attr("transform", "translate(" + main_marginLegend.left + "," + main_marginLegend.top + ")")
		.style("fill", "url(#legendGradient)");

	// append legend text
	main_svg.append("text")
		.attr("text-anchor", "middle")
		.attr("x", -50)
		.attr("y", 15)
		.attr("dy", 0)
		.attr("class", "label legendGradientText")
		.attr("transform", "rotate(-90)")
		.text("Tolerant")
		.style('pointer-events', 'none')
		.style('user-select', 'none');

	// append legend text
	main_svg.append("text")
		.attr("text-anchor", "middle")
		.attr("x", -150)
		.attr("y", 15)
		.attr("dy", 0)
		.attr("class", "label legendGradientText")
		.attr("transform", "rotate(-90)")
		.text("Neutral")
		.style('pointer-events', 'none')
		.style('user-select', 'none');

	// append legend text
	main_svg.append("text")
		.attr("text-anchor", "middle")
		.attr("x", -250)
		.attr("y", 15)
		.attr("dy", 0)
		.attr("class", "label legendGradientText")
		.attr("transform", "rotate(-90)")
		.text("Intolerant")
		.style('pointer-events', 'none')
		.style('user-select', 'none');
}

//Draw the legend for the MetaDomain landscape
function createMetaDomainLegend(){
	// append colors
	main_svg.append("rect")
		.attr("width", 70)
		.attr("height", 20)
		.attr("x", 0)
		.attr("y", 20)
		.attr("dy", 35)
		.attr("class", "label legendMetaDomainRect")
		.style("fill", "green");
	
	// append colors
	main_svg.append("rect")
		.attr("width", 70)
		.attr("height", 20)
		.attr("x", 0)
		.attr("y", 105)
		.attr("dy", 35)
		.attr("class", "label legendMetaDomainRect")
		.style("fill", "red");
	
	// append colors
	main_svg.append("rect")
		.attr("width", 70)
		.attr("height", 20)
		.attr("x", 0)
		.attr("y", 190)
		.attr("dy", 35)
		.attr("class", "label legendMetaDomainRect")
		.style("fill", "black");

	// append legend text
	main_svg.append("text")
		.attr("x", 0)
		.attr("y", 20)
		.attr("dy", 35)
		.attr("class", "label legendMetaDomainText")
		.text("gnomAD")
		.style("font-size", "12px")
		.style('pointer-events', 'none')
		.style('user-select', 'none');
	main_svg.append("text")
		.attr("x", 0)
		.attr("y", 35)
		.attr("dy", 35)
		.attr("class", "label legendMetaDomainText")
		.text("missense in")
		.style("font-size", "12px")
		.style('pointer-events', 'none')
		.style('user-select', 'none');
	main_svg.append("text")
		.attr("x", 0)
		.attr("y", 50)
		.attr("dy", 35)
		.attr("class", "label legendMetaDomainText")
		.text("homologues")
		.style("font-size", "12px")
		.style('pointer-events', 'none')
		.style('user-select', 'none');
	
	// append legend text
	main_svg.append("text")
		.attr("x", 0)
		.attr("y", 105)
		.attr("dy", 35)
		.attr("class", "label legendMetaDomainText")
		.text("ClinVar")
		.style("font-size", "12px")
		.style('pointer-events', 'none')
		.style('user-select', 'none');
	main_svg.append("text")
		.attr("x", 0)
		.attr("y", 120)
		.attr("dy", 35)
		.attr("class", "label legendMetaDomainText")
		.text("missense in")
		.style("font-size", "12px")
		.style('pointer-events', 'none')
		.style('user-select', 'none');
	main_svg.append("text")
		.attr("x", 0)
		.attr("y", 135)
		.attr("dy", 35)
		.attr("class", "label legendMetaDomainText")
		.text("homologues")
		.style("font-size", "12px")
		.style('pointer-events', 'none')
		.style('user-select', 'none');

	// append legend text
	main_svg.append("text")
		.attr("x", 0)
		.attr("y", 190)
		.attr("dy", 35)
		.attr("class", "label legendMetaDomainText")
		.text("no alignment")
		.style("font-size", "12px")
		.style('pointer-events', 'none')
		.style('user-select', 'none');
}

/*******************************************************************************
 * Interactive behaviour functions
 ******************************************************************************/

// Allows toggling between the tolerance landscape and metadomain landscape
function toggleToleranceLandscapeOrMetadomainLandscape(){
	// get the tolerance graph
    var tolerance_graph = d3.select("#tolerance_graph");
    var tolerance_legend = d3.select("#legendGradientRect");
    var tolerance_legend_text = d3.selectAll(".legendGradientText");

    // Get the metadomain graph
    var metadomain_graph = d3.select("#metadomain_graph");
    var meta_legend = d3.selectAll(".legendMetaDomainRect");
    var meta_legend_text = d3.selectAll(".legendMetaDomainText");

    switch($('input[name=landscape_checkbox]:checked', '#checkbox_for_landscape').val()){
        case "metadomain_landscape":
            tolerance_graph.style("opacity", 0);
            tolerance_legend.style("opacity", 0);
            tolerance_legend_text.style("opacity", 0);
            metadomain_graph.style("opacity", 1);
            meta_legend_text.style("opacity", 1);
            meta_legend.style("opacity", 1);
            metadomain_graph_visible = true;
            break;
        case "tolerance_landscape":
            tolerance_graph.style("opacity", 1);
            tolerance_legend.style("opacity", 1);
            tolerance_legend_text.style("opacity", 1);
            metadomain_graph.style("opacity", 0);
            meta_legend_text.style("opacity", 0);
            meta_legend.style("opacity", 0);
            metadomain_graph_visible = false;
            break;
	default:
	    break;
    }
}

function draw_position_schematic_protein(d, element){
	var pathogenic_missense_variant_count = 0;
	if (clinvar_variants_visible){
		// count any pathogenic variants at this position
		if (d.values[0].ClinVar != null) {
			pathogenic_missense_variant_count += d.values[0].ClinVar.length;
		}
	}

	var homologous_pathogenic_missense_variant_count = 0;
	if (homologous_clinvar_variants_visible){
		// count pathogenic variants linked via meta-domain relationships
		if (d.values[0].domains != null){
			meta_domain_ids.forEach(domain_id => {
				if (d.values[0].hasOwnProperty('domains') && d.values[0].domains[domain_id] != null){
					homologous_pathogenic_missense_variant_count = d.values[0].domains[domain_id].pathogenic_missense_variant_count;
				}
			});
		}
	}

	// priortize if selected
    if (d.values[0].selected) {
    	d3.select(element).style("fill-opacity", 0.7);
		return 'green';
    }
    
    // if containing pathogenic variants, display it as red
	if (pathogenic_missense_variant_count > 0){
		d3.select(element).style("fill-opacity", 0.7);
		return 'red';
	}
	
	// if containing pathogenic variants, display it as red
	if (homologous_pathogenic_missense_variant_count > 0){
		d3.select(element).style("fill-opacity", 0.7);
		return 'red';
	}
	
	
	else{
		d3.select(element).style("fill-opacity", 0.2);
		return "grey";
	}
}

function toggleClinvarVariantsInProtein(clinvar_checkbox){
	var focusAxis = d3.select("#tolerance_axis");

	clinvar_variants_visible = clinvar_checkbox.checked;
	
	focusAxis.selectAll(".toleranceAxisTick").style("fill", function(d, i) {
		return draw_position_schematic_protein(d, this);
	});
}

function toggleHomologousClinvarVariantsInProtein(clinvar_checkbox){
	var focusAxis = d3.select("#tolerance_axis");

	homologous_clinvar_variants_visible = clinvar_checkbox.checked;
	
	focusAxis.selectAll(".toleranceAxisTick").style("fill", function(d, i) {
		return draw_position_schematic_protein(d, this);
	});
}


// Rescale the landscape for zooming or brushing purposes
function rescaleLandscape(){
    var focus = d3.select("#tolerance_graph");
    focus.selectAll(".area").attr("d", toleranceArea);
    focus.select(".line").attr("d", toleranceLine);
    focus.select(".axis--x").call(main_xAxis).selectAll("text").remove();
    
    var metadomain_landscape = d3.select("#metadomain_graph");
    metadomain_landscape.selectAll(".normal_missense_variant_count")
	.attr("x", function(d) { return main_x(d.values[0].protein_pos - 0.45); })
	.attr("width", function(d, i) {
	    if (d.values[1].protein_pos != d.values[0].protein_pos){
			return main_x(d.values[1].protein_pos-0.6) - main_x(d.values[0].protein_pos);
	    	    } else {
	    		return main_x(d.values[1].protein_pos +1) - main_x(d.values[0].protein_pos);
	    	    }
		});
    metadomain_landscape.selectAll(".pathogenic_missense_variant_count")
    	.attr("x", function(d) { return main_x(d.values[0].protein_pos+0.05); })
    	.attr("width", function(d, i) {
    	    if (d.values[1].protein_pos != d.values[0].protein_pos){
    			return main_x(d.values[1].protein_pos-0.6) - main_x(d.values[0].protein_pos);
    	    	    } else {
    	    		return main_x(d.values[1].protein_pos +1) - main_x(d.values[0].protein_pos);
    	    	    }
    		});

    metadomain_landscape.selectAll(".not_aligned_position_plot")
    	.attr("x", function(d) { return main_x(d.values[0].protein_pos-0.5); })
    	.attr("width", function(d, i) {
    	    if (d.values[1].protein_pos != d.values[0].protein_pos){
    	    	return main_x(d.values[1].protein_pos) - main_x(d.values[0].protein_pos);
	    	    } else {
	    		return main_x(d.values[1].protein_pos +1) - main_x(d.values[0].protein_pos);
	    	    }
    		});

    var focusAxis = d3.select("#tolerance_axis");

    focusAxis.selectAll(".toleranceAxisTick").attr("x", function(d, i) {
	return main_x(d.values[0].protein_pos - 0.5);
    	})
    	.attr("width", function(d, i) {
    	    if (d.values[1].protein_pos != d.values[0].protein_pos){
		return main_x(d.values[1].protein_pos) - main_x(d.values[0].protein_pos);
    	    } else {
    		return main_x(d.values[1].protein_pos +1) - main_x(d.values[0].protein_pos);
    	    }
	});

    focusAxis.selectAll(".toleranceAxisTickBackground").attr("x", function(d, i) {
	return main_x(d.values[0].protein_pos - 0.5);
    	})
    	.attr("width", function(d, i) {
    	    if (d.values[1].protein_pos != d.values[0].protein_pos){
		return main_x(d.values[1].protein_pos) - main_x(d.values[0].protein_pos);
    	    } else {
    		return main_x(d.values[1].protein_pos +1) - main_x(d.values[0].protein_pos);
    	    }
	});
    
    focusAxis.selectAll(".toleranceAxisTickLabel")
    	.attr("x", function(d, i) {
    	    return main_x(d.values[0].protein_pos);
	})
	.attr("text-anchor", "middle")
	.style("opacity", function(d, i) {
		var textwidth = 13;
	    var rectwidth = d3.select('#toleranceAxisRect_' + d.values[0].protein_pos).node().width.animVal.value;
	    if ((textwidth *0.75) >= rectwidth) {
		return 0;
	    } else if ((textwidth * 0.9) >= rectwidth) {
		return 0.1;
	    } else if ((textwidth * 1.0) >= rectwidth) {
		return 0.4;
	    } else if ((textwidth * 1.25) >= rectwidth) {
		return 0.5;
	    } else if ((textwidth * 1.5) >= rectwidth) {
		return 0.7;
	    }
	    return 1;
	});

    var domains = d3.select("#domain_annotation");
    domains.select(".axis--x").call(main_xAxis);
    domains.selectAll(".pfamDomains")
    	.attr("x", function(d) {
    	    return main_x(d.start - 0.5);
    	})
    	.attr("width", function(d) {
    	    return main_x(d.stop + 1) - main_x(d.start);
	});
    domains.selectAll(".clinvar")
    	.attr("x1", function(d) { 
    	    return main_x(d.protein_pos);
	})
	.attr("x2", function(d) {
	    return main_x(d.protein_pos);
	});
}

//This function resets the zooming
function resetZoom(){
	d3.select("#brush_for_zooming").call(brush.move, null);
}

/*******************************************************************************
 * Static interface functions
 ******************************************************************************/

// the color coding for specific tolerance scores
// color #f29e2e indicates the average dn/ds tolerance score over all genes
function tolerance_color(score) {
	if (score <= 0.175) {
		return toleranceColorGradient[0].color;
	} else if (score <= 0.35) {
		return toleranceColorGradient[1].color;
	} else if (score <= 0.525) {
		return toleranceColorGradient[2].color;
	} else if (score <= 0.7) {
		return toleranceColorGradient[3].color;
	} else if (score <= 0.875) {
		return toleranceColorGradient[4].color;
	} else if (score <= 1.025) {
		return toleranceColorGradient[5].color;
	} else if (score <= 1.2) {
		return toleranceColorGradient[6].color;
	} else if (score <= 1.375) {
		return toleranceColorGradient[7].color;
	} else {
		return toleranceColorGradient[8].color;
	}
}

//the color coding for specific tolerance scores
//color #f29e2e indicates the average dn/ds tolerance score over all genes
function tolerance_rating(score) {
	if (score <= 0.175) {
		return '<label style="background-color:'+tolerance_color(score)+';color:black">highly intolerant</label>';
	} else if (score <= 0.35) {
		return '<label style="background-color:'+tolerance_color(score)+';color:black">intolerant</label>' ;
	} else if (score <= 0.525) {
		return '<label style="background-color:'+tolerance_color(score)+';color:black">intolerant</label>' ;
	} else if (score <= 0.7) {
		return '<label style="background-color:'+tolerance_color(score)+';color:black">slightly intolerant</label>' ;
	} else if (score <= 0.875) {
		return '<label style="background-color:'+tolerance_color(score)+';color:black">neutral</label>' ;
	} else if (score <= 1.025) {
		return '<label style="background-color:'+tolerance_color(score)+';color:black">slightly tolerant</label>' ;
	} else if (score <= 1.2) {
		return '<label style="background-color:'+tolerance_color(score)+';color:black">tolerant</label>' ;
	} else if (score <= 1.375) {
		return '<label style="background-color:'+tolerance_color(score)+';color:black">tolerant</label>' ;
	} else {
		return '<label style="background-color:'+tolerance_color(score)+';color:black">highly tolerant</label>' ;
	}
}