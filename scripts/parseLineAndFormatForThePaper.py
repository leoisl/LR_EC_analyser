def formatValues(line, goodCondition, averageCondition, suffix="\\%", decimalPrecision=2, divideBy=1, toInt=False):
	"""
	Example:
	roundValues("13.723	1.855	1.317	0.439	4.503	3.731	1.592	5.451	4.357	0.381	0.675	2.655	0.33	6.427	5.199	4.122	2.912	4.488	5.647	5.71", 1, 3)
	"""
	values = [round(float(x)/float(divideBy), decimalPrecision) for x in line.split()]
	if toInt:
		for i,value in enumerate(values):
			values[i] = int(value)
	toPrint = []
	for value in values:
		if eval(goodCondition.format(value=value)):
			toPrint.append("\\good{%s%s}"%(value, suffix))
		elif eval(averageCondition.format(value=value)):
			toPrint.append("\\average{%s%s}"%(value, suffix))
		else:
			toPrint.append("\\bad{%s%s}"%(value, suffix))

	print " & ".join(toPrint)




"""
Some command-lines I used in the paper
nb reads:
formatValues("740776        740776  708629  914491  740776  676616  1388145 619139  619139  619172  1321299 738224  626272", "{value}>=641 and {value}<=841", "{value}>=541 and {value}<=941",suffix="k", decimalPrecision=0, divideBy=1000, toInt=True)

mapped reads:
formatValues("83.473	88.108	95.64	98.817	85.512	95.47	97.489	97.099	97.567	98.712	99.204	85.526	98.914", "{value}>=98", "{value}>=92",suffix="\\%", decimalPrecision=1)

mean length:
formatValues("2010.96174282	2173.79600986	1925.68633657	1377.94825318	2096.76706454	1953.09001413	815.632734333	2212.39870368	1900.7684268	1930.57283275	775.509438061	2117.17601433	1796.29380205", "{value}>=1800 and {value}<=2200", "{value}>=1400 and {value}<=2600",suffix="", decimalPrecision=0, toInt=True)

nb bases:
formatValues("1312707425	1468545289	1333662059	1245224564	1393812378	1288587985	1106348407	1331724664	1151046074	1178676189	1015240493	1400254687	1112405437", "{value}>=1213 and {value}<=1413", "{value}>=1000 and {value}<=1600", suffix="M", decimalPrecision=0, divideBy=1000000, toInt=True)

mapped bases:
formatValues("89.018	90.295	96.551	99.16	90.632	95.937	99.122	90.945	97.66	97.517	99.235	92.372	99.486", "{value}>=98", "{value}>=94",suffix="\\%", decimalPrecision=1)

nb detected genes:
nb reads:
formatValues("16816	17037	16818	16550	16772	16581	16478	16530	16186	15040	15572	16580	14632", "{value}>=16 and {value}<=18", "{value}>=13 and {value}<=20",suffix="k", decimalPrecision=1, divideBy=1000)


"""