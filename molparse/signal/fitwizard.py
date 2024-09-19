import dash
import numpy as np
import pandas as pd
import plotly.graph_objects as go
from dash import dcc
from dash import html

from scipy.optimize import curve_fit

from dash.dependencies import Input, Output, State

import plotly.graph_objects as go

from IPython.display import display

# selected_points = []

_logfile = None

# from dash import Dash, dcc, html, Input, Output

# app = Dash(__name__)
app = dash.Dash(prevent_initial_callbacks=True)


# @app.callback(
# 	Output(component_id='my-output', component_property='children'),
# 	Input(component_id='graph', component_property='value')
# )
# def update_output_div(input_value):
# 	return f'Output: {input_value}'


def linear(x, a, b):
    return a * x + b


def fitwizard(xdata, ydata, type="linear", mode="markers", logfile=None):
    global _logfile

    _logfile = logfile

    if _logfile is not None:
        f = open(_logfile, "w")
        f.write("# type x_min x_max slope intercept\n")
        f.close()

    # create the plotly graph

    fig = go.FigureWidget(
        [
            go.Scatter(x=xdata, y=ydata, mode=mode),
            go.Scatter(x=None, y=None, mode="lines"),
            go.Scatter(x=None, y=None, mode="markers"),
        ]
    )

    scatter = fig.data[0]

    # colors = ['#a3a7e4'] * 1

    # scatter.marker.color = colors

    # scatter.marker.size = [10] * 100

    fig.layout.hovermode = "closest"

    # fig.update_layout(
    #    updatemenus=[
    #        dict(
    #            type="buttons",
    #            direction="right",
    #            active=0,
    #            x=0.57,
    #            y=1.2,
    #            buttons=list([
    #                dict(label="Select",
    #                     method="update",
    #                     args=[{"dragmode": "select"},
    #                     ]
    #                     ),
    #                dict(label="Zoom",
    #                     method="update",
    #                     args=[{"dragmode": "zoom"},
    #                     ]
    #                     )
    #                ]
    #                )
    #            )
    #        ]
    #        )

    # fig.update_layout(dragmode="select")

    app.layout = html.Div(
        [
            # html.H6("Change the value in the text box to see callbacks in action!"),
            # html.Div([
            dcc.Graph(id="graph", figure=fig),
            # ]),
            # html.Br(),
            # html.Div(id='my-output'),
        ]
    )

    # color = np.zeros(len(ydata), dtype='uint8')
    # colorscale = [[0, '#167b7e'], [1, '#4b3268']]

    # app = dash.Dash(prevent_initial_callbacks=True)
    # app.layout = html.Div([dcc.Graph(figure=build_figure(), id="graph")])

    app.run_server(debug=True)


@app.callback(
    Output("graph", "figure"),
    [Input("graph", "selectedData"), Input("graph", "clickData")],
    [State("graph", "figure")],
)
def linear_fit(selectedData, clickData, fig):
    global _logfile

    selection = None
    # Update selection based on which event triggered the update.
    trigger = dash.callback_context.triggered[0]["prop_id"]
    if trigger == "graph.clickData":
        selection = [point["pointNumber"] for point in clickData["points"]]
    if trigger == "graph.selectedData":
        selection = [point["pointIndex"] for point in selectedData["points"]]

    # Update scatter selection
    fig["data"][0]["selectedpoints"] = selection

    xdata = np.array(fig["data"][0]["x"])
    ydata = np.array(fig["data"][0]["y"])

    sel_x = xdata[selection]
    sel_y = ydata[selection]

    param, covariance = curve_fit(linear, sel_x, sel_y)

    fig["data"][1] = go.Scatter(
        x=sel_x, y=linear(sel_x, param[0], param[1]), mode="lines"
    )

    print(f"Fitting in xrange [{min(sel_x)},{max(sel_x)}]")

    print(f"Slope: {param[0]}")
    print(f"Intercept: {param[1]}")

    if _logfile is not None:
        f = open(_logfile, "a")
        f.write(f"LIN {min(sel_x)} {max(sel_x)} {param[0]} {param[1]}\n")
        f.close()

    return fig
