"""
    PlotTool: Your Go-To Python Plotting Package
    ---
        since: 230920
        auth: Andy6

    ---
    `PlotTool` combines popular `Python plotting tools` into one package.
    ---
    Goal: 
    ---
        Make plotting easy and fast.

    ---
    How? :
    ---
        By using consistent parameter names.

    ---
    Who is it for?
    ---
        Anyone who wants quick and simple visualizations.
	---
    quick start:
    ---
        basic arguments:
        ---
        - data : input data (DataFrame)
        - x, y : column names of data (str)
            - For plot bar: x and y axis
        - x_label, y_label : label names of plot x and y axis (str)
        - title : plot title (str)
        - percentage : Whether to display the percentage of each bar relative to the total (bool)
"""

from . import plt
