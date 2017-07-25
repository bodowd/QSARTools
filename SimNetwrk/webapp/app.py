from flask import Flask, render_template
from flask_bootstrap import Bootstrap

app = Flask(__name__)
Bootstrap(app)
@app.route('/mol_network/')
def mol_network():
    return render_template('mol_network.html')
if __name__ == '__main__':
    app.run(debug=True)
