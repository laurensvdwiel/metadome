/**
 * Controls all responsiveness for the dashboard.html page
 */

// Initialize with global_debugging set to false
var global_debugging = false;

function setGenomeBuildValue(genomeBuild) {
    const genomeBuildElement = document.getElementById("genomeBuild");
    if (!genomeBuildElement) {
        return;
    }

    genomeBuildElement.value = genomeBuild;
    genomeBuildElement.innerText = genomeBuild;
    genomeBuildElement.textContent = genomeBuild;
}

function getGenomeBuildValue() {
    const genomeBuildElement = document.getElementById("genomeBuild");
    if (!genomeBuildElement) {
        return "";
    }

    return genomeBuildElement.value || genomeBuildElement.innerText || genomeBuildElement.textContent || "";
}

function setGeneNameValue(geneName) {
    const geneNameElement = document.getElementById("geneName");
    if (!geneNameElement) {
        return;
    }
    geneNameElement.value = geneName || "";
}

// Specify the style for the console messages
const info_style = "background: blue; color: white; padding-left:10px;";
const success_style = "background: green; color: white; padding-left:10px;";

// This function allows clicking events to be raised on d3js elements
jQuery.fn.d3Click = function () {
  this.each(function (i, e) {
    var evt = new MouseEvent("click");
    e.dispatchEvent(evt);
  });
};

$("#geneName").keyup(function(event) {
	// Allow the user to press enter to retrieve transcripts
    if (event.keyCode === 13) {
        $("#getTranscriptsButton").click();
    }
});

// when the user presses the ESC the overlays are disabled
document.addEventListener("keydown", keyPress, false);
function keyPress (e) {
    if (global_debugging) {console.info('%c Function keyPress(e) called with e = ' + e, info_style);}
    if(e.key === "Escape") {
    	$("#domain_information_overlay").removeClass("is-active");
    	$("#positional_information_overlay").removeClass("is-active");
        d3.selectAll(".tr").classed("is-selected", false);

        // New: if the URL ends with /p.<n>/, remove that segment (shareable state cleared)
        try {
            const url = new URL(window.location.href);
            const segments = url.pathname.split('/').filter(s => s.length > 0);
            if (segments.length > 0 && /^p\.\d+$/.test(segments[segments.length - 1])) {
                segments.pop();
                url.pathname = '/' + segments.join('/') + '/';
                window.history.replaceState({}, "", url.toString());
            }
        } catch (err) {
            // no-op
        }
    }
}

function get_query_param(param) {
    if (global_debugging) {console.info('%c Function get_query_paran(param) called with param = ' + param, info_style);}
	var result =  window.location.search.match(
			new RegExp("(\\?|&)" + param + "(\\[\\])?=([^&]*)"));
    return result ? result[3] : false;
}

function show_tour_data_input(index) {
    if (global_debugging) {console.info('%c Function show_tour_data_input(index) called with index = ' + index, info_style);}
    switch (index) {
		case 0:// reset
			tour_reset_gene_name();
			break;
		case 1: // genome build
			break;
		case 2: // gene input
			// ensure no interaction is possible
			document.getElementById("retrieve_transcripts_control").setAttribute("style", "pointer-events: none;");
			document.getElementById("start_analysis_control").setAttribute("style", "pointer-events: none;");

			//reset
			document.getElementById("geneNameHelpMessage").innerHTML = "";

			// behaviour
			tour_fill_gene_name();
			break;
		case 3: // messages after get transcripts
			//reset
			resetDropdown();

			// behaviour
			tour_fill_succes_message_transcripts();
			break;
		case 4: // start analysis
			// reset
			$("#loading_overlay").removeClass('is-active');

			// behaviour
			tour_fill_transcripts();
			break;
		case 5: // loading overlay
			// reset
			$("#toleranceGraphContainer").addClass('is-hidden');
			$("#graph_control_field").addClass('is-hidden');
			$("#selected_positions_information").addClass('is-hidden');

			// behaviour
			$("#loading_overlay").addClass('is-active');
			$("#loading_label").text("Loading...");
			break;
		case 6: // fill the graph
			$("#loading_overlay").removeClass('is-active');
			tour_fill_graph_example();
			break;
		case 8:
			// resets
			resetZoom();
			break;
		case 9: // zooming explanation
			// resets
			resetZoom();
			tour_turn_clinvar_variants_off();
			tour_check_default_visualization();toggleToleranceLandscapeOrMetadomainLandscape();
			document.getElementById("graph_control_field_checkboxes").setAttribute("style", "pointer-events: none;");
			break;
		case 10: // landscape modes
			// resets
			resetZoom();
			tour_turn_clinvar_variants_off();

			// behaviour
			window.setTimeout(tour_check_alternative_visualization, 500);
			window.setTimeout(toggleToleranceLandscapeOrMetadomainLandscape,1000);
			window.setTimeout(tour_check_default_visualization, 2200);
			window.setTimeout(toggleToleranceLandscapeOrMetadomainLandscape,2200);
			document.getElementById("graph_control_field_checkboxes").removeAttribute("style", "pointer-events: none;");
			break;
		case 11: // clinvar variants
			// resets
			document.getElementById("graph_control_field_checkboxes").removeAttribute("style", "pointer-events: none;");
			resetZoom();

			// behaviour
			tour_switch_clinvar_variants();
			window.setTimeout(tour_switch_clinvar_variants,1200);
			break;
		case 12: // explanation of schematic protein
			// resets
			tour_check_default_visualization();
			tour_turn_clinvar_variants_off();
			tour_unselect_position_in_domain();
			resetZoom();
			document.getElementById("graph_control_field_checkboxes").setAttribute("style", "pointer-events: none;");
			break;
		case 13: // positional information
			// resets
			resetZoom();
			tour_turn_clinvar_variants_off();

			// Behaviour
			window.setTimeout(tour_select_position_in_domain, 500);
			break;
		}
}

function tour_turn_clinvar_variants_off(){
    if (global_debugging) {console.info('%c Function tour_turn_clinvar_variants_off() called', info_style);}
	if ($('#clinvar_checkbox').is(':checked') == true ){
		tour_switch_clinvar_variants();
	}
}

function tour_switch_clinvar_variants(){
    if (global_debugging) {console.info('%c Function tour_switch_clinvar_variants() called', info_style);}
	$('#clinvar_checkbox').trigger( "click" );
}

function tour_unselect_position_in_domain(){
    if (global_debugging) {console.info('%c Function tour_unselect_position_in_domain() called', info_style);}
	var tour_selected_position;
	d3.select('#toleranceAxisRect_68').attr("element_is_selected",function(d,i){ tour_selected_position = d.values[0].selected});
	
	if (tour_selected_position){
		tour_trigger_position_of_interest();
	}
}

function tour_select_position_in_domain(){
    if (global_debugging) {console.info('%c Function tour_select_position_in_domain() called', info_style);}
	var tour_selected_position;
	d3.select('#toleranceAxisRect_68').attr("element_is_selected",function(d,i){ tour_selected_position = d.values[0].selected});
	
	if (!tour_selected_position){
		tour_trigger_position_of_interest();
	}
}

function tour_trigger_position_of_interest(){
    if (global_debugging) {console.info('%c Function tour_trigger_position_of_interest() called', info_style);}
	$('#toleranceAxisRect_68').d3Click();
}

function tour_check_default_visualization(){
    if (global_debugging) {console.info('%c Function tour_check_default_visualization() called', info_style);}
	$('#checkbox_for_landscape_default').prop('checked', true);
}
function tour_check_alternative_visualization(){
    if (global_debugging) {console.info('%c Function tour_check_alternative_visualization() called', info_style);}
	$('#checkbox_for_landscape_alternative').prop('checked', true);
}

function tour_default_genome_build(){
    if (global_debugging) {console.info('%c Function tour_default_genome_build() called', info_style);}
	const genome_build_tour = 'GRCh38'
	var genome_build_button = document.getElementById("genomeBuild");
	genome_build_button.innerText = genome_build_tour;
}

function tour_fill_gene_name(){
    if (global_debugging) {console.info('%c Function tour_fill_gene_name() called', info_style);}
	$('#geneName').val('T');
}

function tour_reset_gene_name(){
    if (global_debugging) {console.info('%c Function tour_reset_gene_name() called', info_style);}
	$('#geneName').val('');
}

function tour_fill_succes_message_transcripts(){
    if (global_debugging) {console.info('%c Function tour_fill_succes_message_transcripts() called', info_style);}
	document.getElementById("geneNameHelpMessage").innerHTML = "Retrieved transcripts for gene 'T'";
	$("#geneNameHelpMessage").removeClass('is-danger');
	$("#geneNameHelpMessage").addClass('is-success');
}

function tour_fill_transcripts(){
    if (global_debugging) {console.info('%c Function tour_fill_transcripts() called', info_style);}
	clearDropdown();
	$("#getToleranceButton").addClass('is-info');
	$("#getToleranceButton").removeClass('is-static');
	var dropdown = document.getElementById("gtID");
	dropdown.setAttribute('class', 'dropdown');
	var opt = new Option();
	opt.value = 1;
	opt.text = "ENST00000296946.2 / NM_003181.3 (435aa)";
	dropdown.options.add(opt);
}

function tour_fill_graph_example(){
    if (global_debugging) {console.info('%c Function tour_fill_graph_example() called', info_style);}
	// Retrieve the example json to fill in the graph
	$.getJSON(window.METADOME_CONFIG.exampleGeneJsonUrl, function(json) {
		document.getElementById("toleranceGraphContainer").setAttribute("style", "pointer-events: none;");
		$("#toleranceGraphContainer").removeClass('is-hidden');
	    $("#graph_control_field").removeClass('is-hidden');
	    var geneName = document.getElementById("geneName").value;
	    var geneDetails = document.getElementById("geneDetails");
		geneDetails.innerHTML = 'Gene: '+json.gene_name+' (transcript: <a href="'+ window.METADOME_LINKS.ensemblTranscript(json.genome_build || getGenomeBuildValue(), json.transcript_id) +'" target="_blank">'+json.transcript_id+'</a>, protein: <a href="'+ window.METADOME_LINKS.uniprot(json.protein_ac) +'" target="_blank">'+json.protein_ac+'</a>)';
		createGraph(json);
	    // disable any mouse events on the elements
		document.getElementById("graph_control_field_buttons").setAttribute("style", "pointer-events: none;");
		document.getElementById("graph_control_field_checkboxes").setAttribute("style", "pointer-events: none;");
		document.getElementById("selected_positions_information").setAttribute("style", "pointer-events: none;");
	});
}

//Setup tour
var tour = new Tour({
  backdrop: true,
  orphan: true,
  storage: false,
  onStart: function (tour) {tour_default_genome_build();},
  onNext: function (tour) {
    show_tour_data_input(tour.getCurrentStep() + 1)
  },
  onPrev: function (tour) {
    show_tour_data_input(tour.getCurrentStep() - 1)
  },
  onEnd: function (tour) {  window.location.reload(true); window.location.replace(window.METADOME_CONFIG.dashboardUrl);},
  steps: [
	  {
		element: "",
		title: "Start Tour",
		content: "Welcome to MetaDome.<br><br>"+
				"This tour explains the usage of MetaDome "+
				"via an example.<br><br>"+
				"You can click 'end tour' any time to start "+
				"analysing your own gene of interest. " +
				"<br><br>Click next to begin."
	  },
	  {
	    element: "#genomeBuild",
	    title: "Genome Build",
	    content: "MetaDome support the human genome builds GRCh37 "+
	    		"and GRCh38. By default the newest ('GRCh38') is" +
	    		"selected. You can always check which genome build" +
	    		"is selected in this field."
	  },
	  {
	    element: "#retrieve_transcripts_control",
	    title: "Gene of interest",
	    content: "Input here your gene name of interest Then you" +
	    		" can click the 'Get Transcripts' button to " +
	    		"retrieve all the transcripts for your gene of " +
	    		"interest.<br><br> For this example we fill in " +
	    		"the 'T' gene."
	  },
	  {
		element: "#geneNameHelpMessage",
		title: "Get Transcripts",			
		content: "After clicking the 'Get Transcripts' button, we "+
				"check if there are any transcripts available for "+
				"your gene of interest.<br><br>"+
				"If the gene was not present in our database you " +
				"would have seen the following message:<br><br>" +
				"<div class='help is-danger'>No transcripts " +
				"available in database for gene 'T'</div>"				
	  },
	  {
		element: "#start_analysis_control",
		title: "Select transcript & analyse your protein",
		content: "Select the transcript of your interest from this "+
				"dropdown menu. You can than click the 'Analyse " +
				"Protein' button to start the analysis."
	  },
	  {
		element: "#loading_overlay",
		title: "Loading screen",
		content: "This analysis may take between 5 minutes and up " +
				"to an hour to complete, but in this example: <br><br>" +
				"<b>you can click next to continue</b>. <br><br> " +
				"Luckily all results are stored after they "+
				"complete, so the next time you query the same " +
				"transcript it will load in a matter of seconds."
	  },
	  {
		element: "#content",
		title: "Analysis results",
		content: "When the analysis completes, you obtain a "+
				"wealth of information. <br><br> But don't worry, " +
				"we will now go over each part in detail.",
	  },
	  {
		element: "#graph_control_field",
		title: "Graph Control",
		content: "Here, you may switch between different " +
				"representations of the analysis result visualizations," +
				" download the current view as an '.svg', reset any " +
				"zooming of the graph or select, or display any known " +
				"ClinVar variants.",
	  },
	  {
		element: "#landscape_svg",
		title: "Visualization",
		content: "Here you may find the visualization of your results." +
				"It is a graphic representation of the protein of your " +
				"transcript of interest.<br><br> By default the " +
				"metadomain variants are displayed.",
	  },
	  {
		element: "#schematic_protein_zoom_text",
		title: "Zooming in",
		content: "You can zoom in on regions of your interest by " +
				"clicking anywhere in the gray area and dragging " +
				"the mouse. If you just click on any part in the " +
				"protein the zooming is reset.<br><br> Go ahead " +
				"and try for yourself",
		placement: 'left',
		backdropContainer: "#landscape_svg"
	  },
	  {
		element: "#checkbox_for_landscape",
		title: "Switching between visualizations",
		content: "Here you can switch between the 'Meta-domain " +
				"Landscape' and the 'Tolerance Landscape' " +
				"visualization modes.<br><br> Go ahead " +
				"and try for yourself",
		backdropContainer: "#landscape_svg"
	  },
	  {
		element: "#clinvar_checkbox",
		title: "Switching ClinVar variants",
		content: "Here you can toggle the display of any ClinVar " +
				"variants for your gene of interest. <br><br> Go ahead " +
				"and try for yourself",
		backdropContainer: "#landscape_svg"
	  },
	  {
		element: "#tolerance_axis",
		title: "Schematic protein representation",
		content: "For each gene a schematic protein is displayed for " +
				"all the positions and the presence of any Pfam protein " +
				"domains is annotated.<br><br> Here each position is " +
				"hoverable and selectable. If you click a position you " +
				"obtain much more information about it.",
		backdropContainer: "#landscape_svg"
	  },
	  {
		element: "#toleranceAxisRect_68",
		title: "Selecting positions of interest",
		content: "If you click a position it will become highlighted. " +
				"And you may view more detailed information for that " +
				"position.",
		backdropContainer: "#landscape_svg"
	  },
	  {
		element: "#selected_positions_information",
		title: "More info for a selected position",
		content: "All the highlighted positions that you have selected " +
				"will be put into this list. Here a short overview is " +
				"displayed on any known gnomAD or ClinVar variants " +
				"present at this position. Also variants that are " +
				"homologously related to this position are displayed " +
				"if the position is part of a Pfam protein domain." +
				"<br><br> Clicking a selected position provides a " +
				"pop-up with even more details.",
	  },
]});
tour.init();

// Function to handle receiving all transcripts belonging to a gene name from the webserver
function getTranscript(geneName, genomeBuild, transcript_id_results) {
    if (global_debugging) {console.info('%c Function getTranscript(geneName, genomeBuild, transcript_id_results) called, with geneName='+geneName+', genomeBuild='+genomeBuild+', transcript_id_results='+ transcript_id_results, info_style);}

    // Clear any previous transcripts in user interface
    clearTranscripts();

    // Check if there are any transcripts retrieved
    if (typeof geneName !== 'undefined' && geneName.length > 0) {

        // Check if there are any messages to be displayed
        document.getElementById("geneNameHelpMessage").innerHTML = transcript_id_results.message;

        if (transcript_id_results.transcript_ids.length > 0) {
            clearDropdown();
            $("#geneName").addClass('is-success');
            $("#geneName").removeClass('is-danger');
            $("#geneNameHelpMessage").addClass('is-success');
            $("#geneNameHelpMessage").removeClass('is-danger');
            $("#getToleranceButton").addClass('is-info');
            $("#getToleranceButton").removeClass('is-static');
            var dropdown = document.getElementById("gtID");
            dropdown.setAttribute('class', 'dropdown');

            // Sort the results by sequence length (longest first)
            transcript_id_results.transcript_ids.sort(function(a,b){return (a.aa_length<b.aa_length) ? 1 : ((b.aa_length<a.aa_length) ? -1: 0);});

            var maneSelectIndex = -1; // Track MANE_Select transcript index

            // Add the transcripts as options
            for (var i = 0; i < transcript_id_results.transcript_ids.length; i++) {
                var opt = new Option();
                opt.value = i;

                // Build the display text
                var displayText = transcript_id_results.transcript_ids[i].gencode_id;

                if (transcript_id_results.transcript_ids[i].refseq_nm_numbers !== "") {
                    displayText += " / " + transcript_id_results.transcript_ids[i].refseq_nm_numbers;
                }

                displayText += " (" + transcript_id_results.transcript_ids[i].aa_length + "aa)";

                // Add MANE transcript type if available and not empty
                if (transcript_id_results.transcript_ids[i].mane_transcript_type &&
                    transcript_id_results.transcript_ids[i].mane_transcript_type !== "") {
                    displayText += " [" + transcript_id_results.transcript_ids[i].mane_transcript_type + "]";

                    // Check if this is MANE_Select (note the capital S)
                    if (transcript_id_results.transcript_ids[i].mane_transcript_type === "MANE_Select") {
                        maneSelectIndex = i;
                    }
                }

                opt.text = displayText;

                if (!transcript_id_results.transcript_ids[i].has_protein_data) {
                    opt.disabled = true;
                }

                dropdown.options.add(opt);
            }

            // Select MANE_Select transcript by default if available and has protein data
            if (maneSelectIndex !== -1 && transcript_id_results.transcript_ids[maneSelectIndex].has_protein_data) {
                dropdown.selectedIndex = maneSelectIndex;
            } else {
                // Fallback to first available transcript with protein data
                for (var i = 0; i < transcript_id_results.transcript_ids.length; i++) {
                    if (transcript_id_results.transcript_ids[i].has_protein_data) {
                        dropdown.selectedIndex = i;
                        break;
                    }
                }
            }

        } else {
            $("#geneName").addClass('is-danger');
            $("#geneName").removeClass('is-success');
            $("#geneNameHelpMessage").addClass('is-danger');
            $("#geneNameHelpMessage").removeClass('is-success');
            $("#getToleranceButton").addClass('is-static');
            $("#getToleranceButton").removeClass('is-info');
            resetDropdown();
        }
    }
}

// function to clear all items in the transcript dropdown
function clearDropdown() {
    if (global_debugging) {console.info('%c Function clearDropdown() called', info_style);}
	var dropdown = document.getElementById("gtID");
	for (var i = dropdown.options.length - 1; i >= 0; i--) {
		dropdown.remove(i);
	}
	dropdown.disabled = false;
}

function resetDropdown() {
    if (global_debugging) {console.info('%c Function resetDropdown() called', info_style);}
	clearDropdown()

	var dropdown = document.getElementById("gtID");
	dropdown.disabled = true;

	var opt = new Option();
	opt.value = "";
	opt.disabled = true;
	opt.selected = true;
	opt.hidden = true;
	opt.text = "Please first retrieve the transcripts";
	dropdown.options.add(opt);
}

function resetGraphControl(){
    if (global_debugging) {console.info('%c Function resetGraphControl() called', info_style);}
	document.getElementById("clinvar_checkbox").checked = false;
	document.getElementById("checkbox_for_landscape_default").checked = true;
}

function clearTranscripts() {
    if (global_debugging) {console.info('%c Function clearTranscripts() called', info_style);}
	resetDropdown();
	resetGraph();
	resetGraphControl();
	$("#getToleranceButton").addClass('is-static');
	$("#getToleranceButton").removeClass('is-info');
	$("#geneName").removeClass('is-success');
	$("#advanced_options").addClass("is-hidden");
	$("#toleranceGraphContainer").addClass('is-hidden');
    	document.getElementById("geneNameHelpMessage").innerHTML = "";
	$("#graph_control_field").addClass('is-hidden');
}

function saveSvg(svgEl, name) {
    if (global_debugging) {console.info('%c Function saveSvg(svgEl, name) called', info_style);}
    svgEl.setAttribute("xmlns", "http://www.w3.org/2000/svg");
    var svgData = svgEl.outerHTML;
    var preface = '<?xml version="1.0" standalone="no"?>\r\n';
    var svgBlob = new Blob([preface, svgData], {type:"image/svg+xml;charset=utf-8"});
    var svgUrl = URL.createObjectURL(svgBlob);
    var downloadLink = document.createElement("a");
    downloadLink.href = svgUrl;
    downloadLink.download = name+'.svg';
    document.body.appendChild(downloadLink);
    downloadLink.click();
    document.body.removeChild(downloadLink);
}

function visualize() {
	if (global_debugging) {console.info('%c Function visualize() called', info_style);}
    resetGraphControl();
    var selection = document.getElementsByClassName("dropdown")[0];

    if (selection !== undefined && selection !== null && selection.length > 0) {
        var input = selection.options[selection.selectedIndex].text;
        var gtID = input.toUpperCase();

        if (gtID !== undefined && gtID.length > 0 && !selection.options[selection.selectedIndex].disabled) {
            var transcript_id = gtID.split(" ")[0];
            var genome_build = getGenomeBuildValue();

            if (!genome_build) {
                console.error("visualize(): genome build is empty");
                return;
            }

            $("#loading_overlay").addClass('is-active');
            $("#loading_label").text("Loading...");
            $.ajax(
                {
                    type: 'POST',
                    url: window.METADOME_CONFIG.submitVisualizationUrl,
                    data: JSON.stringify({'transcript_id': transcript_id, 'genome_build': genome_build}),
                    success: function(data) { getVisualizationStatus(data.transcript_id, data.genome_build, 0); },
                    error: function(response) {
                        obj = JSON.parse(response.responseText);
                        $("#error-feedback").text(obj['error']);
                        $("#loading_overlay").removeClass('is-active');
                    },
                    contentType: "application/json",
                    dataType: 'json'
                }
            );
        }
    }
}

function getVisualizationStatus(transcript_id, genome_build, checkCount) {
    if (global_debugging) {console.info('%c Function getVisualizationStatus(transcript_id, checkCount) called, with transcript_id='+transcript_id+', genome_build='+genome_build+', checkCount='+checkCount, info_style);}
    $.get(
        Flask.url_for("api.get_visualization_status_for_transcript",
            {'transcript_id': transcript_id, 'genome_build': genome_build}
        ),
        {
            _: Date.now()
        },
        function(data) {
            if (data.status == 'SUCCESS')
                getVisualizationResult(transcript_id, genome_build);
            else if (data.status == 'FAILURE')
                showVisualizationError(transcript_id, genome_build);
            else {   // try again after a while..
                var delay = 10000;
                if (checkCount >= 5) {
                    delay = 50000;
                    $("#loading_label").html("Loading...</br>This is taking longer than usual, and may take up to an hour. You can choose to wait or check back later.");
                }
                checkCount++;

                setTimeout(function() { getVisualizationStatus(transcript_id, genome_build, checkCount); }, delay);
            }
        }
    );
}

function getVisualizationResult(transcript_id, genome_build) {
	if (global_debugging) {
		console.info('%c Function getVisualizationResult(transcript_id) called, with transcript_id=' + transcript_id + ' genome_build=' + genome_build, info_style);
	}
	$.get(
		Flask.url_for("api.get_visualization_result_for_transcript",
			{'transcript_id': transcript_id, 'genome_build': genome_build}
		),
		{
			_: Date.now()
		},
		function (obj) {
			$("#loading_overlay").removeClass('is-active');
			$("#toleranceGraphContainer").removeClass('is-hidden');
			$("#graph_control_field").removeClass('is-hidden');
			$("#geneName").val(obj.gene_name);

			var geneDetails = document.getElementById("geneDetails");

			var refSeqLinks = "";
			const refseqIdsRaw = obj.refseq_ids;
			const refseqArray = Array.isArray(refseqIdsRaw)
				? refseqIdsRaw.map(id => String(id || "").trim()).filter(id => id.length > 0)
				: String(refseqIdsRaw || "").split(',').map(id => id.trim()).filter(id => id.length > 0);

			if (refseqArray.length > 0) {
				refSeqLinks += "RefSeq: ";

				for (var i = 0; i < refseqArray.length; i++) {
					if (i > 0) {
						refSeqLinks += ', ';
					}
					refSeqLinks += '<a href="' + window.METADOME_LINKS.refseq(refseqArray[i]) + '" target="_blank">' + refseqArray[i] + '</a>';
				}
				refSeqLinks += ',';
			}

			geneDetails.innerHTML = 'Protein of ' + obj.gene_name + ' (GENCODE: <a href="'
				+ window.METADOME_LINKS.ensemblTranscript(genome_build, obj.transcript_id) + '" target="_blank">' + obj.transcript_id
				+ '</a>, ' + refSeqLinks + ' UniProt: <a href="'
				+ window.METADOME_LINKS.uniprot(obj.protein_ac) + '" target="_blank">' + obj.protein_ac + '</a>)';
			createGraph(obj);

			d3.select('#dlSVG').on('click', function () {
				var fileName = 'Gene_' + obj.gene_name + '_transcript_' + obj.transcript_id + '_protein_' + obj.protein_ac + '_metadome';
				saveSvg(document.getElementById('landscape_svg'), fileName);
			});

			checkForSelectedPosition();
		}
	);
}

function showVisualizationError(transcript_id, genome_build) {
    if (global_debugging) {console.info('%c Function showVisualizationError(transcript_id) called, with transcript_id='+transcript_id+' genome_build='+genome_build, info_style);}
    // Forward to the error page.
    customRedirect(routeUrl("web.visualization_error", {'transcript_id': transcript_id, 'genome_build': genome_build}));
}

// annotates meta domain information for a position
function createPositionalInformation(domain_metadomain_coverage, transcript_id, position_json) {
    if (global_debugging) {console.info('%c Function createPositionalInformation(omain_metadomain_coverage, transcript_id, position_json) called', info_style);}
    // Retrieve the needed information for the GET request
    var protein_position = position_json.values[0].protein_pos;

    // Construct the request for this domain and the aligned positions
    var requested_domains = {};
    var domain_ids = Object.keys(position_json.values[0].domains);
    for (var i = 0; i < domain_ids.length; i++){
        var domain_id = domain_ids[i];

        // Check if this position is meta domain suitable
        if (!(position_json.values[0].domains[domain_id] == null)){		
            requested_domains[domain_id] = [];

            // append the consensus positions
            for (var j = 0; j < position_json.values[0].domains[domain_id].consensus_pos.length; j ++) {
                requested_domains[domain_id].push(position_json.values[0].domains[domain_id].consensus_pos[j])
            }
        }
    }

    if (Object.keys(requested_domains).length > 0) {
        // Activate the loading overlay
        $("#loading_overlay").addClass('is-active');
        $("#loading_label").text("Loading...");

        var genome_build = document.getElementById("genomeBuild").value;

        var input = {"requested_domains": requested_domains,
                     "transcript_id": transcript_id,
                     "protein_position": protein_position,
                     "genome_build": genome_build};

        // Execute the POST request
        $.ajax(
          {
             type: 'POST',
             url: window.METADOME_CONFIG.metadomainAnnotationUrl,
             data: JSON.stringify(input),
             success:function(data) {
                $("#loading_overlay").removeClass('is-active');
                FillPositionalInformation(domain_metadomain_coverage, position_json, data);
                $("#positional_information_overlay").addClass('is-active');
             },
             contentType: "application/json",
             dataType: 'json'
          }
        );
    }
    else {
        // No domains requested, so we can fill in the information without performing a GET request
        FillPositionalInformation(domain_metadomain_coverage, position_json, {});
        $("#positional_information_overlay").addClass('is-active');
    }
}

//Adds positional information for a selected position
function FillPositionalInformation(domain_metadomain_coverage, position_data, data){
    if (global_debugging) {console.info('%c Function FillPositionalInformation(domain_metadomain_coverage, position_data, data) called', info_style);}
	// Reset the positional information
    document.getElementById("positional_information_overlay_title").innerHTML = '<div class="label"><label class="title">Positional information (p.'+ position_data.values[0].protein_pos+')</label></div>';
    document.getElementById("positional_information_overlay_body").innerHTML = '';

    
    // Add information on position to the HTML
    document.getElementById("positional_information_overlay_body").innerHTML += '<label class="label">Protein details</label>';
    document.getElementById("positional_information_overlay_body").innerHTML += '<p>'+document.getElementById("geneDetails").innerHTML +'</p>';
    document.getElementById("positional_information_overlay_body").innerHTML += '</br><label class="label">Location details</label>';
    document.getElementById("positional_information_overlay_body").innerHTML += '<p>Chr: '+position_data.values[0].chr+', strand: '+position_data.values[0].strand+'</p>';
    document.getElementById("positional_information_overlay_body").innerHTML += '<p>Gene: '+ position_data.values[0].chr_positions +'</p>';
    document.getElementById("positional_information_overlay_body").innerHTML += '<p>Protein: p.'+ position_data.values[0].protein_pos +' '+ position_data.values[0].ref_aa_triplet+'</p>';
    document.getElementById("positional_information_overlay_body").innerHTML += '<p>cDNA: '+ position_data.values[0].cdna_pos +' '+ position_data.values[0].ref_codon +'</p>';
    document.getElementById("positional_information_overlay_body").innerHTML += '<p>Exon: '+ position_data.values[0].exon_numbers +'</p>';
    document.getElementById("positional_information_overlay_body").innerHTML += '<p>Tolerance score (dn/ds): '+ (Math.round((position_data.values[0].sw_dn_ds)*100)/100) +' ('+tolerance_rating(position_data.values[0].sw_dn_ds) +')</p>';
    
    
    // retrieve domain information
    var domain_information = "";
    var meta_domain_information = "";
    if (Object.keys(position_data.values[0].domains).length > 0){
		var domain_ids = '';
		var domain_id_list = Object.keys(position_data.values[0].domains);
		var n_domains_at_position = Object.keys(position_data.values[0].domains).length;
		
		for (var i = 0; i < n_domains_at_position; i++){
		    if (i+1 == n_domains_at_position){
			domain_ids += '<a href="' + window.METADOME_LINKS.pfam(domain_id_list[i]) + '" target="_blank">'+domain_id_list[i]+"</a>";
		    }
		    else{
			domain_ids += '<a href="' + window.METADOME_LINKS.pfam(domain_id_list[i]) + '" target="_blank">'+domain_id_list[i]+", </a>";
		    }
		    
		    // add meta domain information
		    if (!(position_data.values[0].domains[domain_id_list[i]] == null)){
		    	// retrieve the domain_information from the data object
		    	var meta_domain = data[domain_id_list[i]];
		    	
		    	// compute coverage
		    	var position_coverage = Math.round(((meta_domain.alignment_depth/domain_metadomain_coverage[domain_id_list[i]])*100)*10)/10;
		    	
		    	// Add information to the report
				meta_domain_information += '<hr><label class="label">Meta-domain information for domain '+domain_id_list[i]+':</label>';
				meta_domain_information += '<p>Aligned to consensus position '+ position_data.values[0].domains[domain_id_list[i]].consensus_pos+', related to '+ (meta_domain.alignment_depth-1) +' other codons throughout the genome (with a '+position_coverage+'\% alignment coverage).</p>';
				
				
				var gnomAD_table = '<hr><label class="label">Variants in gnomAD SNVs at homologous positions:</label>';
				gnomAD_table += createGnomADTableHeader();
				var clinvar_table = '<hr><label class="label">Pathogenic ClinVar SNVs at homologous positions:</label>';
				clinvar_table += createClinVarTableHeader();
				
				// Append gnomad
				gnomAD_table += createGnomADTableBody(meta_domain.normal_variants);
				
				// Append clinvar
		    	clinvar_table += createClinVarTableBody(meta_domain.pathogenic_variants);
				
				// Add the footers
				clinvar_table += createTableFooter();
				gnomAD_table += createTableFooter();
				
				// Reset the tables if there are no variants
				if (meta_domain.normal_variants.length == 0){
				    gnomAD_table = "";
				}
				if (meta_domain.pathogenic_variants.length == 0){
				    clinvar_table = "";
				}
				
				// Add the meta-domain information to the domain information
				meta_domain_information += clinvar_table + gnomAD_table;
		    }
		    else{
			meta_domain_information += '<label class="label">No meta-domain information at position for domain '+domain_id_list[i]+'</label>';
		    }
		}
		
		// Add the domain info to the html element
		domain_information += '<p>Position is part of protein domain(s): '+domain_ids+'</p>'
    }
    
    // Add the domain information to the html element
    document.getElementById("positional_information_overlay_body").innerHTML += domain_information

    // Add clinvar at position information
    document.getElementById("positional_information_overlay_body").innerHTML += '<hr>';
	document.getElementById("positional_information_overlay_body").innerHTML += '<label class="label">Known pathogenic ClinVar SNVs at position</label>';
    if ("ClinVar" in position_data.values[0]){
		// Add ClinVar variant table
    	var clinvar_variants = position_data.values[0].ClinVar;
    	
    	// retrieve the gene name
        var geneName = document.getElementById("geneName").value;
    	
    	// Add the positional information
        for (var index = 0; index < clinvar_variants.length; index++){        	
        	clinvar_variants[index].gene_name = geneName;
        	clinvar_variants[index].chr = position_data.values[0].chr;
        	clinvar_variants[index].chr_positions = position_data.values[0].chr_positions;
        	clinvar_variants[index].ref_codon = position_data.values[0].ref_codon;
        	clinvar_variants[index].ref_aa_triplet = position_data.values[0].ref_aa_triplet;
        }
        
        // Add the clinvar variant information
    	document.getElementById("positional_information_overlay_body").innerHTML += createClinVarTableHeader()+ createClinVarTableBody(clinvar_variants)+ createTableFooter();
    }
    else{
    	document.getElementById("positional_information_overlay_body").innerHTML += '<p>No ClinVar SNVs found at position</p>';
    }
    
    // Add the meta-domain information to the html element
    document.getElementById("positional_information_overlay_body").innerHTML += meta_domain_information;
}

//Update the positional information table with new values
function addRowToPositionalInformationTable(domain_metadomain_coverage, d, transcript_id) {
    if (global_debugging) {console.info('%c Function addRowToPositionalInformationTable(domain_metadomain_coverage, d, transcript_id) called', info_style);}
	var new_row = d3.select('#position_information_tbody').append('tr').attr('class', 'tr').attr("id", "positional_table_info_" + d.values[0].protein_pos);
	
	new_row.append('th').text(d.values[0].protein_pos);
	new_row.append('td').text(d.values[0].ref_aa_triplet);

	var domain_ids = "-";
	var clinvar_at_pos = "-";
	var related_gnomad = "-";
	var related_clinvar = "-";
	
	// Add clinvar at position information
	if ("ClinVar" in d.values[0]){
	    clinvar_at_pos = ""+d.values[0].ClinVar.length;
	}
	else{
	    clinvar_at_pos = "0";
	}
	
	// add domain and metadomain information to the information
	if (Object.keys(d.values[0].domains).length > 0){
	    var domain_id_list = Object.keys(d.values[0].domains);
	    var n_domains_at_position = Object.keys(d.values[0].domains).length;
	    domain_ids = "";
	    related_gnomad = 0;
	    related_clinvar = 0;
	    for (var i = 0; i < n_domains_at_position; i++){
			if (i+1 == n_domains_at_position){
			    domain_ids += domain_id_list[i];
			}	
			else{
			    domain_ids += domain_id_list[i]+", ";
			}
			// append normal and pathogenic variant count
			if (!(d.values[0].domains[domain_id_list[i]] == null)){
	    		    related_gnomad += d.values[0].domains[domain_id_list[i]].normal_variant_count;
	    		    related_clinvar += d.values[0].domains[domain_id_list[i]].pathogenic_variant_count;
			}
			else{
				related_gnomad = "-";
				related_clinvar = "-";
			}
	    }
	}
	new_row.append('td').text(domain_ids);
	new_row.append('td').text(clinvar_at_pos);
	new_row.append('td').text(related_gnomad);
	new_row.append('td').text(related_clinvar);
	
	// Add interactiveness to the rows
	new_row.on("click", function() {
	    d3.selectAll('.tr').classed("is-selected", false);
	    d3.select(this).classed("is-selected", true);
	    
	    // Call this method found in dashboard.js
	    createPositionalInformation(domain_metadomain_coverage, transcript_id, d)
	}).on("mouseover", function(d, i) {
	    d3.select(this).style("cursor", "pointer");
	});
		
	// Sort the table to the protein positions
	sortTable();
}

function createClinVarTableHeader(){
    if (global_debugging) {console.info('%c Function createClinVarTableHeader() called', info_style);}
    var html_table= '';
    // Define the header
    html_table += '<table class="table is-hoverable is-narrow">';
    html_table += '<thead><tr style="border-style:hidden;">';
    html_table += '<th><abbr title="Gene Symbol">Gene</abbr></th>';
    html_table += '<th><abbr title="Chromosomal position">Position</abbr></th>';
    html_table += '<th><abbr title="Change of nucleotide">Variant</abbr></th>';
    html_table += '<th><abbr title="Change of residue">Residue change</abbr></th>';
    html_table += '<th><abbr title="Type of mutation">Type</abbr></th>';
    html_table += '<th><abbr title="ClinVar Identifier">ClinVar ID</abbr></th>';
    html_table += '</tr></thead><tfoot></tfoot><tbody>';
    return html_table;
}

function createClinVarTableBody(ClinvarVariants){
    if (global_debugging) {console.info('%c Function createClinVarTableBody(ClinvarVariants) called', info_style);}
    var html_table= '';
    // here comes the data
    for (var index = 0; index < ClinvarVariants.length; index++){
		var variant = ClinvarVariants[index];
		html_table += '<tr>';
		html_table += '<td>'+variant.gene_name+'</td>';
		html_table += '<td><a href="' + window.METADOME_LINKS.ensemblGenomicPosition(variant.genome_build || getGenomeBuildValue(), variant.chr, variant.pos) + '" target="_blank">'+variant.chr+':'+variant.pos+'</a></td>';
		html_table += '<td>'+variant.ref+'>'+variant.alt+'</td>';
		html_table += '<td>'+variant.ref_aa_triplet+'>'+variant.alt_aa_triplet+'</td>';
		html_table += '<td>'+variant.type+'</td>';
		html_table += '<td><a href="' + window.METADOME_LINKS.clinvarVariant(variant.clinvar_ID) + '" target="_blank">' + variant.clinvar_ID + '</a></td>';
		html_table += '</tr>';
    }
    	
    return html_table;
}

function createGnomADTableHeader(){
    if (global_debugging) {console.info('%c Function createGnomADTableHeader() called', info_style);}
    var html_table= '';
    // Define the header
    html_table += '<table class="table is-hoverable is-narrow">';
    html_table += '<thead><tr style="border-style:hidden;">';
    html_table += '<th><abbr title="Gene Symbol">Gene</abbr></th>';
    html_table += '<th><abbr title="Chromosomal position">Position</abbr></th>';
    html_table += '<th><abbr title="Change of nucleotide">Variant</abbr></th>';
    html_table += '<th><abbr title="Change of residue">Residue change</abbr></th>';
    html_table += '<th><abbr title="Type of mutation">Type</abbr></th>';
    html_table += '<th><abbr title="Allele frequency">gnomAD Allele Frequency</abbr></th>';
    html_table += '</tr></thead><tfoot></tfoot><tbody>';
    return html_table;
}

function createGnomADTableBody(gnomADVariants){
    if (global_debugging) {console.info('%c Function createGnomADTableBody(gnomADVariants) called', info_style);}
    var html_table= '';
    // here comes the data
    for (var index = 0; index < gnomADVariants.length; index++){
		var variant = gnomADVariants[index];
		html_table += '<tr>';
		html_table += '<td>'+variant.gene_name+'</td>';
		html_table += '<td><a href="' + window.METADOME_LINKS.ensemblGenomicPosition(variant.genome_build || getGenomeBuildValue(), variant.chr, variant.pos) + '" target="_blank">'+variant.chr+':'+variant.pos+'</a></td>';
		html_table += '<td>'+variant.ref+'>'+variant.alt+'</td>';
		html_table += '<td>'+variant.ref_aa_triplet+'>'+variant.alt_aa_triplet+'</td>';
		html_table += '<td>'+variant.type+'</td>';
		html_table += '<td><a href="' + window.METADOME_LINKS.gnomadVariant(variant.genome_build || getGenomeBuildValue(), variant.chr, variant.pos, variant.ref, variant.alt) + '" target="_blank">'+parseFloat(variant.allele_count/variant.allele_number).toFixed(6)+'</a></td>';
		html_table += '</tr>';
    }
    	
    return html_table;
}

function createTableFooter(){
    if (global_debugging) {console.info('%c Function createTableFooter() called', info_style);}
    return '</tbody></table>';
}

function sortTable() {
    if (global_debugging) {console.info('%c Function sortTable() called', info_style);}
	var rows = $('#position_information_tbody tr').get();
	
	rows.sort(function(a, b) {
	    var A = parseInt($(a).children('th').eq(0).text());
	    var B = parseInt($(b).children('th').eq(0).text());
	
	    if(A < B) {
	      return -1;
	    }
	
	    if(A > B) {
	      return 1;
	    }
	
	    return 0;
	});
	
	$.each(rows, function(index, row) {
	  $('#position_information_tbody').append(row);
    });
}

function customRedirect(url){
    // redirect while keeping debug settings
    if (global_debugging) {
		window.location.assign(url + "?debug=True");
	} else {
		window.location.assign(url);
	}
}

function routeUrl(endpoint, params) {
    if (window.Flask && typeof Flask.url_for === "function") {
        return Flask.url_for(endpoint, params || {});
    }

    // fallback using config
    if (endpoint === "web.dashboard") return window.METADOME_CONFIG.dashboardUrl;
    if (endpoint === "web.dashboard_genome_build") {
        return window.METADOME_CONFIG.dashboardBaseUrl + "/" + encodeURIComponent(params.genome_build) + "/";
    }
    if (endpoint === "web.dashboard_gene_name") {
        return window.METADOME_CONFIG.dashboardBaseUrl + "/" + encodeURIComponent(params.genome_build) + "/" + encodeURIComponent(params.gene_name) + "/";
    }
    if (endpoint === "web.transcript") {
        return window.METADOME_CONFIG.dashboardBaseUrl + "/" + encodeURIComponent(params.genome_build) + "/" + encodeURIComponent(params.gene_name) + "/" + encodeURIComponent(params.transcript_id) + "/";
    }

    throw new Error("Unknown endpoint fallback: " + endpoint);
}

function getTranscriptButtonPress() {
    if (global_debugging) {console.info('%c Function getTranscriptButtonPress() called', info_style);}
		const geneName = document.getElementById("geneName").value;
		const genomeBuild = getGenomeBuildValue();

		if (typeof genomeBuild !== 'undefined' && genomeBuild.length > 0) {
			if (typeof geneName !== 'undefined' && geneName.length > 0) {
				customRedirect(routeUrl("web.dashboard_gene_name", {'genome_build': genomeBuild, 'gene_name': geneName}));
			}
			else {
				customRedirect(routeUrl("web.dashboard_genome_build", {'genome_build': genomeBuild}));
			}
		}
		else {
            console.info('%c No genome build specified in the URL, redirecting to the default dashboard.', info_style);
            customRedirect(routeUrl("web.dashboard"))
		}
}

function getToleranceButtonPress() {
    if (global_debugging) {console.info('%c Function getToleranceButtonPress() called', info_style);}
	const geneName = document.getElementById("geneName").value;
	const genomeBuild = getGenomeBuildValue();
	const selection = document.getElementsByClassName("dropdown")[0];

	if (typeof genomeBuild !== 'undefined' && genomeBuild.length > 0) {
		if (typeof geneName !== 'undefined' && geneName.length > 0) {
			if (selection !== undefined && selection !== null && selection.length > 0) {
				var input = selection.options[selection.selectedIndex].text;
				var gtID = input.toUpperCase();

				if (gtID !== undefined && gtID.length > 0 && !selection.options[selection.selectedIndex].disabled) {
					const transcript_id = gtID.split(" ")[0];
					customRedirect(routeUrl("web.transcript", {'genome_build': genomeBuild, 'gene_name': geneName, 'transcript_id': transcript_id}));
				}
			}
		}
	}
}

async function checkForSelectedPosition() {
    if (global_debugging) {console.info('%c Function checkForSelectedPosition() called', info_style);}

    const path = window.location.pathname;
    const pathSegments = path.split('/').filter(segment => segment.length > 0);

    if (pathSegments.length < 6) {
        if (global_debugging) {
            console.info('%c No position segment present in path', info_style);
        }
        return;
    }

    const positionSegment = pathSegments[5];
    const match = /^p\.(\d+)$/.exec(positionSegment);

    if (!match) {
        if (global_debugging) {
            console.info('%c No valid p.<position> segment found in path', info_style);
        }
        return;
    }

    const position = parseInt(match[1], 10);
    const elementId = 'toleranceAxisRect_' + position;

    let attempts = 0;
    const maxAttempts = 100;

    const checkElement = setInterval(function() {
        attempts++;
        const element = document.getElementById(elementId);

        if (element) {
            clearInterval(checkElement);

            const transcript_id = element.getAttribute('data-transcript-id');
            const domain_coverage_str = element.getAttribute('data-domain-coverage');

            if (!transcript_id || !domain_coverage_str) {
                console.error('Missing data attributes on element');
                return;
            }

            const domain_metadomain_coverage = JSON.parse(domain_coverage_str);
            const d3Element = d3.select('#' + elementId);
            const boundData = d3Element.datum();

            if (boundData && transcript_id && domain_metadomain_coverage) {
                if (!boundData.values[0].selected) {
                    d3Element.style("fill", "green").style("fill-opacity", 0.7);
                    boundData.values[0].selected = true;

                    addRowToPositionalInformationTable(
                        domain_metadomain_coverage,
                        boundData,
                        transcript_id
                    );

                    selected_positions += 1;
                    $("#position_information_table").removeClass('is-hidden');
                    document.getElementById("selected_positions_explanation").innerHTML =
                        'Click on one of the selected positions in the table to view more information';
                }

                createPositionalInformation(
                    domain_metadomain_coverage,
                    transcript_id,
                    boundData
                );
            }
        } else if (attempts >= maxAttempts) {
            clearInterval(checkElement);
            if (global_debugging) {
                console.info('%c Position p.' + position + ' not found after ' + maxAttempts + ' attempts', info_style);
            }
        }
    }, 100);
}

function parseDashboardPath() {
    if (global_debugging) {console.info('%c Function parseDashboardPath() called', info_style);}

    const path = window.location.pathname;
    let pathSegments = path.split('/').filter(segment => segment.length > 0);

    if (get_query_param('debug')) {
        global_debugging = true;
        console.info('%c Debugging enabled', success_style);
        console.info('%c Path: ' + path, info_style);
        console.info('%c Path segments length: ' + pathSegments.length, info_style);
        console.info('%c Path segments: ' + pathSegments, info_style);
    }

    if (!(pathSegments.length > 1 && pathSegments[0] === 'metadome' && pathSegments[1] === 'dashboard')) {
        console.info('%c Invalid URL path: "' + path + '", expected format is /metadome/dashboard/<GenomeBuild>/<gene_name>/<transcript_id>.', info_style);
        customRedirect(routeUrl("web.dashboard"));
        return null;
    }

    const parsed = {
        genome_build: pathSegments[2] || "",
        gene_name: pathSegments[3] || "",
        transcript_id: pathSegments[4] || "",
        position: null,
        is_tour: pathSegments[2] === 'tour'
    };

    if (pathSegments.length > 5) {
        const match = /^p\.(\d+)$/.exec(pathSegments[5]);
        if (match) {
            parsed.position = parseInt(match[1], 10);
        }
    }

    if (global_debugging) {
        console.info('%c Parsed dashboard path:', info_style);
        console.info(parsed);
    }

    return parsed;
}

async function fetchTranscriptResults(geneName, genomeBuild) {
    if (global_debugging) {
        console.info('%c Function fetchTranscriptResults(geneName, genomeBuild) called with geneName=' + geneName + ', genomeBuild=' + genomeBuild, info_style);
    }

    try {
        const response = await fetch(
            Flask.url_for("api.get_transcript_ids_for_gene", {
                'genome_build': genomeBuild,
                'gene_name': geneName
            }),
            {
                cache: "no-store"
            }
        );

        if (global_debugging) {
            console.info('%c Transcript fetch response status: ' + response.status, info_style);
        }

        if (!response.ok) {
            throw new Error('get_transcript_ids_for_gene did not respond: ' + response.statusText);
        }

        const result = await response.json();

        if (global_debugging) {
            console.info('%c Transcript fetch result:', info_style);
            console.info(result);
        }

        return result;
    } catch (error) {
        console.error('Error fetching transcript data: ' + error);
        return {
            transcript_ids: [],
            message: "No transcripts available in database for gene '" + geneName + "'",
            genome_build: genomeBuild,
            gene_name: geneName
        };
    }
}

async function handleURLPathParameters() {
    if (global_debugging) {console.info('%c Function handleURLPathParameters() called', info_style);}

    const urlState = parseDashboardPath();
    if (urlState === null) {
        console.error("URL state parsing returned null");
        return;
    }

    if (global_debugging) {
        console.info('%c URL state before handling:', info_style);
        console.info(urlState);
    }

    switch (true) {
        case urlState.is_tour:
            if (global_debugging) { console.info('%c Tour URL detected', info_style); }
            tour.restart();
            return;

        case !!urlState.genome_build && !urlState.gene_name:
            if (global_debugging) { console.info('%c Case: genome build only', info_style); }
            setGenomeBuildValue(urlState.genome_build);
            return;

        case !!urlState.genome_build && !!urlState.gene_name && !urlState.transcript_id:
            if (global_debugging) { console.info('%c Case: genome build + gene name', info_style); }

            setGeneNameValue(urlState.gene_name);
            setGenomeBuildValue(urlState.genome_build);

            const transcriptResultsForGene = await fetchTranscriptResults(urlState.gene_name, urlState.genome_build);
            getTranscript(urlState.gene_name, urlState.genome_build, transcriptResultsForGene);

            if (global_debugging) {
                console.info('%c After getTranscript(), dropdown length = ' + document.getElementById("gtID").options.length, info_style);
            }
            return;

        case !!urlState.genome_build && !!urlState.gene_name && !!urlState.transcript_id:
            if (global_debugging) { console.info('%c Case: genome build + gene name + transcript id', info_style); }

            setGeneNameValue(urlState.gene_name);
            setGenomeBuildValue(urlState.genome_build);

            const transcriptResultsForTranscript = await fetchTranscriptResults(urlState.gene_name, urlState.genome_build);
            getTranscript(urlState.gene_name, urlState.genome_build, transcriptResultsForTranscript);

            var dropdown = document.getElementById("gtID");
            var transcriptFound = false;

            if (global_debugging) {
                console.info('%c Dropdown options after getTranscript(): ' + dropdown.options.length, info_style);
            }

            for (var i = 0; i < dropdown.options.length; i++) {
                if (global_debugging) {
                    console.info('%c Checking transcript ID: ' + dropdown.options[i].text + ' against ' + urlState.transcript_id, info_style);
                }
                if (dropdown.options[i].text.startsWith(urlState.transcript_id)) {
                    dropdown.selectedIndex = i;
                    transcriptFound = true;
                    break;
                }
            }

            if (!transcriptFound) {
                console.error("Transcript from URL not found in dropdown: " + urlState.transcript_id);
                return;
            }

            if (dropdown.options[dropdown.selectedIndex].disabled) {
                console.error("Transcript from URL is present but disabled: " + urlState.transcript_id);
                return;
            }

            visualize();
            return;

        default:
            console.info('%c Invalid URL path state, redirecting to dashboard.', info_style);
            customRedirect(routeUrl("web.dashboard"));
    }
}

document.addEventListener("DOMContentLoaded", function() {
    handleURLPathParameters().catch(function(error) {
        console.error(error);
    });
});