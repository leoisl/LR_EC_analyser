import plotly
from plotly.graph_objs import Scatter, Layout

with open("test.html", "w") as htmlFile:
    htmlFile.write('<html><head><script src="https://cdn.plot.ly/plotly-latest.min.js"></script></head><body>%s</body></html>'%plotly.offline.plot([Scatter(x=[1, 2, 3, 4], y=[4, 3, 2, 1])], include_plotlyjs=False, output_type='div'))

