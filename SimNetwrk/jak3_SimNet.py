import flask

app = flask.Flask(__name__, static_folder = 'jak3_SimNet')

@app.route('/<path:path>')
def static_proxy(path):
    return app.send_static_file(path)
print('Go to http://localhost:8000/jak3_SimNet.html')
app.run(port=8000)
