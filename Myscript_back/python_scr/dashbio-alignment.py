import sys
#import six.moves.urllib.request as urlreq
#from six import PY3

import dash
import dash_bio as dashbio
import dash_html_components as html

if len(sys.argv) < 2:
    print("Usage: python3 {} aligned.fasta".format(sys.argv[0]))
    exit()


with open(sys.argv[1], 'r') as fh:
    data = fh.read()

external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']

app = dash.Dash(__name__, external_stylesheets=external_stylesheets)

#if PY3:
#    data = data.decode('utf-8')

app.layout = html.Div([
    dashbio.AlignmentChart(
        id='my-alignment-viewer',
        data=data
    ),
    html.Div(id='alignment-viewer-output')
])


@app.callback(
    dash.dependencies.Output('alignment-viewer-output', 'children'),
    [dash.dependencies.Input('my-alignment-viewer', 'eventDatum')]
)
def update_output(value):
    if value is None:
        return 'No data.'
    return str(value)


if __name__ == '__main__':
    app.run_server(debug=True)
