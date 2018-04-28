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

