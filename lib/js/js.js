//a global variable holding all HOTs
var HOTs = {};


//checks if the fileserver is running or not
function isSiteOnline(url,callback) {
    var timer = setTimeout(function(){
        // timeout after 5 seconds
        callback(false);
    },5000)

    var img = document.createElement("img");
    img.onload = function() {
        clearTimeout(timer);
        callback(true);
    }

    img.onerror = function() {
        clearTimeout(timer);
        callback(false);
    }

    img.src = url+"a.jpg";
}

function createIGVLink(IGVInfo) {
  return "<a href='lib/html/igv.html?"+IGVInfo+"' target=_blank>Open</a>";
}

function IGVRenderer (instance, td, row, col, prop, value, cellProperties) {
    var IGVInfo = Handsontable.helper.stringify(value);
    td.innerHTML = createIGVLink(IGVInfo);
    return td;
}

function unitInCharToDivisor(unitInChar) {
	switch (unitInChar) {
		case 'M':
		case 'm':
			return 1000000
		case 'K':
		case 'k':
			return 1000
		case 'N':
		case 'n':
			return 1
		default:
			throw "Unknown unitInChar: " + unitInChar
	}
}

function getUnitAsChar(valueRadio) {
	if (valueRadio=='N' || valueRadio=='n')
		return ""
	else
		return valueRadio
}

function niceNumberRenderer (instance, td, row, col, prop, value, cellProperties) {
	//check if we are dealing with percentages, so no need to change
	metric = instance.getDataAtCell(row,0)
	if (metric.indexOf('%') > -1) {
		//no need to change
		td.innerHTML = value.toString() + "%"
		return td;
	}

    //dealing with raw numbers
	idRadio = instance.rootElement.id+"_number_format"
	valueRadio = $("input[name=" + idRadio + "]:checked").val()
	divisor = unitInCharToDivisor(valueRadio)
	valueToShow = value/divisor
	td.innerHTML = parseFloat(valueToShow.toFixed(1)).toString() + getUnitAsChar(valueRadio)
	return td;
}

function fillFeatureTable(feature, tableId, rawData, tools) {

        var tableColumns = [
            {type: 'text'}, //feature name
            {type: 'numeric'}, //largest discrepancy
            {renderer: IGVRenderer} //the IGV renderer
        ]
        for (var i = 0; i < tools.length; i++)
          tableColumns.push({type: 'numeric'})

        var tableColHeaders = [
            feature,
            'Delta',
            'See it in IGV'
        ].concat(tools)


    var featureTableSettings = {
        data: rawData,
        columns: tableColumns,
        colHeaders: tableColHeaders,
        //colWidths: [300, 50],
        copyColsLimit: 1000000,
        copyRowsLimitNumber: 1000000,
        readOnly: true,
        stretchH: 'all',
        wordWrap: false,
        allowInsertColumn: false,
        allowInsertRow: false,
        allowRemoveColumn: false,
        allowRemoveRow: false,
        autoColumnSize: {useHeaders: true},
        autoWrapCol: true,
        autoWrapRow: true,
        manualColumnResize: true,
        columnSorting: true,
        sortIndicator: true
    };
    var featureTableContainer = document.getElementById(tableId);
    HOTs[tableId] = new Handsontable(featureTableContainer, featureTableSettings);
    HOTs[tableId].sort(1, false);
}

function fillMainSubtable(tableId, rawData, tools) {

    var tableColumns=[{type: 'text'}]
    for (var i = 0; i < tools.length; i++)
      tableColumns.push({renderer: niceNumberRenderer})

    var tableColHeaders = ['Metric'].concat(tools)


    var tableSettings = {
        data: rawData,
        columns: tableColumns,
        colHeaders: tableColHeaders,
        //colWidths: [300, 50],
        copyColsLimit: 1000000,
        copyRowsLimitNumber: 1000000,
        readOnly: true,
        stretchH: 'all',
        wordWrap: false,
        allowInsertColumn: false,
        allowInsertRow: false,
        allowRemoveColumn: false,
        allowRemoveRow: false,
        autoColumnSize: {useHeaders: true},
        autoWrapCol: true,
        autoWrapRow: true,
        manualColumnResize: true,
        columnSorting: true,
        sortIndicator: true
    };
    var tableContainer = document.getElementById(tableId);
    HOTs[tableId] = new Handsontable(tableContainer, tableSettings);
}


function setUpIGV(serverURL, divId) {
    decodedURI = decodeURIComponent(window.location.search)
    parsedQuery = querystring.parse(decodedURI)

    $(document).ready(function () {
        trackLists=[]
        trackLists.push({
                name: "Genes",
                type: "annotation",
                sourceType: "file",
                url: serverURL+parsedQuery.annotationURL,
                displayMode: "EXPANDED",
                visibilityWindow: parsedQuery.visibilityWindow
        })
        for (var i=0; ; i++) {
            bam="bam_"+i
            tool="tool_"+i
            if (bam in parsedQuery) {
                trackLists.push({
                    url: serverURL+parsedQuery[bam],
                    name: parsedQuery[tool],
                    visibilityWindow: parsedQuery.visibilityWindow
                })
            }else
                break;
        }


        var div = $("#" + divId)[0],
                options = {
                    showNavigation: true,
                    showRuler: true,
                    reference: {
                        fastaURL: serverURL+parsedQuery.fastaURL
                    },
                    locus: parsedQuery.locus,
                    tracks: trackLists
                };

        igv.createBrowser(div, options);

    });
}


function showDataOnModal(name, data) {
    if ( plotInfo[name] !== null) {
        xLabel = data.points[0].x
        indexTool = data.points[0].curveNumber
        randomString = Math.random().toString(36).substring(2, 15) + Math.random().toString(36).substring(2, 15)
        $('<div id="'+randomString+'"></div>').appendTo('body');
        $('#'+randomString).html("<textarea readonly rows='40' style='width:100%; font-family: \"Courier New\", Courier, \"Lucida Sans Typewriter\", \"Lucida Typewriter\", monospace; font-size: 16px'>"+plotInfo[name][xLabel][indexTool]+"</textarea>")
        $('#'+randomString).iziModal({onClosed: function(){$('#'+randomString).remove()}});
        $('#'+randomString).iziModal('open');
    }
}

function blockUI() {
    $.blockUI({
      message: '<img width="25px" src="lib/resources/busy.gif" /> Loading, please wait...' ,
      css: {
            border: 'none',
            padding: '15px',
            backgroundColor: '#000',
            '-webkit-border-radius': '10px',
            '-moz-border-radius': '10px',
            opacity: .5,
            color: '#fff'
        }
      });
}

function unblockUI() {
    $.unblockUI({
        fadeOut: 0
    })
}