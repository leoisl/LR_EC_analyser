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
  return "<a href='igv.html?"+IGVInfo+"' target=_blank>Open</a>";
}

function IGVRenderer (instance, td, row, col, prop, value, cellProperties) {
    var IGVInfo = Handsontable.helper.stringify(value);
    td.innerHTML = createIGVLink(IGVInfo);
    return td;
}

function fillGeneTable(rawData, tools) {

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
        console.log(tableColumns)

        var tableColHeaders = [
            'Gene',
            'Coord',
            'Strand'
        ].concat(tools).concat(['Largest discrepancy', 'See it in IGV'])


    var geneTableSettings = {
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
    var geneTableContainer = document.getElementById('LR_EC_gene_table');
    geneTable = new Handsontable(geneTableContainer, geneTableSettings);
    geneTable.sort(tableColHeaders.length-1, false); //change this
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
                order: Number.MAX_VALUE,
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
