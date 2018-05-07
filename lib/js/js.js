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
	idRadio = instance.rootElement.id+"_number_format"
	valueRadio = $("input[name=" + idRadio + "]:checked").val()
	divisor = unitInCharToDivisor(valueRadio)
	valueToShow = value/divisor
	td.innerHTML = parseFloat(valueToShow.toFixed(2)).toString() + getUnitAsChar(valueRadio)
    return td;
}

function fillFeatureTable(feature, tableId, rawData, tools) {

        var table3FirstTypeColumns = [
            {type: 'text'},
            {type: 'text'},
            {type: 'text'}
        ]
        var tableToolsTypeColumns=[]
        for (var i = 0; i < tools.length; i++)
          tableToolsTypeColumns.push({type: 'numeric'})
        tableToolsTypeColumns.push({type: 'numeric'}) //largest discrepancy
        tableToolsTypeColumns.push({renderer: IGVRenderer}) //the IGV renderer
        var tableColumns = table3FirstTypeColumns.concat(tableToolsTypeColumns)

        var tableColHeaders = [
            feature,
            'Coord',
            'Strand'
        ].concat(tools).concat(['Largest discrepancy', 'See it in IGV'])


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
    HOTs[tableId].sort(tableColHeaders.length-2, false);
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
                //TODO
                //indexURL: "https://s3.amazonaws.com/igv.broadinstitute.org/annotations/hg19/genes/refGene.hg19.bed.gz.tbi",
                indexed: false, //TODO: check this
                visibilityWindow: 300000000,
                displayMode: "EXPANDED"
        })
        for (var i=0; ; i++) {
            bam="bam_"+i
            tool="tool_"+i
            if (bam in parsedQuery) {
                trackLists.push({
                    url: serverURL+parsedQuery[bam],
                    name: parsedQuery[tool]
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
