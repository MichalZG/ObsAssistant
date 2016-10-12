<!DOCTYPE html>
<html lang="en">
  <head>

    <meta charset="utf-8">
    <meta name="viewport" content="width=device-width, initial-scale=1, shrink-to-fit=no">
    <meta http-equiv="x-ua-compatible" content="ie=edge">

      	<link rel="stylesheet" type="text/css" href="http://dc-js.github.io/dc.js/css/dc.css"/>
        <script type="text/javascript" src="http://dc-js.github.io/dc.js/js/d3.js"></script>
        <script type="text/javascript" src="http://dc-js.github.io/dc.js/js/crossfilter.js"></script>
        <script type="text/javascript" src="http://dc-js.github.io/dc.js/js/dc.js"></script>

    	<link rel="stylesheet" type="text/css" href="css/bootstrap.min.css">
    	<link rel="stylesheet" type="text/css" href="css/style.css">
</head>

<body>

<!-- ******* HEADER SECTION ********* -->
<header>
	<div class="container">

	<h2>ObsAssistant <small>1.0.0</small></h2>

	</div>

</header>

<!-- ******* MAIN SECTION ********* -->
<div class="charts-section">
	<div class="container">

		<div class="row">

			<h3>Charts</h3>

			<div id="chart-ring-year"></div>
			<div id="chart-row-spenders"></div>

		</div>

	</div>
</div>

<!-- ******* FOOTER SECTION ********* -->
<footer>
	<div class="container">
			<div class="row">

				<div class="col-lg-6 col-md-6">

					<p>Program for show in real time during the observations FLUX, SNR, FWHM and position of observed object. Program can help centered object on image (in future this will be automatically), make good focus, watch good SNR etc.</p>

				</div>

				<div class="col-lg-6 col-md-6">

					Footer 2

				</div>

			</div>

			<div class="row">
				<div class="copy">
				Copyright 2016
				</div>
		</div>
	</div>
</footer>


<script type="text/javascript">
    var yearRingChart   = dc.pieChart("#chart-ring-year"),
        spenderRowChart = dc.rowChart("#chart-row-spenders");
    var connection = new WebSocket('ws://localhost:8001/websocket');
    var data1 = [
        {Name: 'Ben', Spent: 330, Year: 2014, 'total':1},
        {Name: 'Aziz', Spent: 1350, Year: 2012, 'total':2},
        {Name: 'Vijay', Spent: 440, Year: 2014, 'total':2},
        {Name: 'Jarrod', Spent: 555, Year: 2015, 'total':1},
    ];
    // set crossfilter with first dataset
    var xfilter = crossfilter(data1),
        yearDim  = xfilter.dimension(function(d) {return +d.Year;}),
        spendDim = xfilter.dimension(function(d) {return Math.floor(d.Spent/10);}),
        nameDim  = xfilter.dimension(function(d) {return d.Name;}),
      
        spendPerYear = yearDim.group().reduceSum(function(d) {return +d.Spent;}),
        spendPerName = nameDim.group().reduceSum(function(d) {return +d.Spent;});
    function render_plots(){
        yearRingChart
            .width(200).height(200)
            .dimension(yearDim)
            .group(spendPerYear)
            .innerRadius(50);
        spenderRowChart
            .width(250).height(200)
            .dimension(nameDim)
            .group(spendPerName)
            .elasticX(true);
        dc.renderAll();
    }
    render_plots();
    // data reset function (adapted)
    function resetData(ndx, dimensions) {
        var yearChartFilters = yearRingChart.filters();
        var spenderChartFilters = spenderRowChart.filters();
        yearRingChart.filter(null);
        spenderRowChart.filter(null);
        xfilter.remove();
        yearRingChart.filter([yearChartFilters]);
        spenderRowChart.filter([spenderChartFilters]);
    }
    connection.onmessage = function(event) {
        var newData = JSON.parse(event.data);
        var updateObject =[{
            "Name": newData.Name,
            "Year": newData.Year,
            "Spent": newData.Spent,
            "payType": newData.payType
        }]
        //resetData(ndx, [yearDim, spendDim, nameDim]);
        xfilter.add(updateObject);
        dc.redrawAll();
    }
</script>

</body>
</html>