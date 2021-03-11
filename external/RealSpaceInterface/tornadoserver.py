from Calc2D.CalculationClass import Calculation

import time
import numpy as np
from concurrent.futures import ThreadPoolExecutor
from tornado.ioloop import IOLoop
from tornado import gen
import tornado.web
import tornado.websocket
import os
import os.path
import json
import unicodedata
import logging
import base64
import traceback
import sys

import config

pool = ThreadPoolExecutor(max_workers=config.MAX_THREADPOOL_WORKERS)

def generate_redshifts(redshift_config):
    logging.info(redshift_config)
    arrs = []
    for conf in redshift_config:
        log = conf["log"]
        func = np.logspace if log else np.linspace
        start = np.log10(conf["from"]) if log else conf["from"]
        stop = np.log10(conf["to"]) if log else conf["to"]
        arrs.append(func(start, stop, conf["points"]))
    # Remove duplicates
    return np.flip(np.unique(np.concatenate(arrs)), axis=0)

# Load available colormaps
def get_colormaps(path=config.COLORMAP_PATH):
    categories = []
    maps = []
    order = {'Default': 1, 'Uniform': 2, 'Diverging': 3, 'Miscellaneous': 4}
    cmap_directories = list(sorted(
        os.listdir(os.path.join("static", path)),
        key=lambda d: order[d]
        ))
    for directory in cmap_directories:
        categories.append(directory)
        maps_for_category = []
        for cmap in os.listdir(os.path.join("static", path, directory)):
            maps_for_category.append({
                'label': cmap[:cmap.rfind(".")],
                'src': os.path.join(os.path.join(config.COLORMAP_PATH, directory, cmap)),
                })
        maps.append(maps_for_category)
    return categories, maps

class SimulationHandler(tornado.web.RequestHandler):
    def get(self):
        categories, colormaps = get_colormaps()
        self.render('RSI.html', categories=categories, colormaps=colormaps)

class DataConnection(tornado.websocket.WebSocketHandler):
    def open(self):
        logging.info("Client connected!")
        self.calc = Calculation(kbins=config.TRANSFER_FUNCTION_CLIENT_SAMPLES)
        # Send list of `k` values only once
        logging.info("Sending k range to client");
        self.write_message(json.dumps({
            "type": "krange",
            "k": self.calc.krange.tolist()
            }))

    def on_close(self):
        logging.info("Connection was closed")

    @gen.coroutine
    def on_message(self, message):
        message = json.loads(message)
        param_type = message['type']
        logging.debug("Received message from client: {}".format(message))
        params = message['params']
        if param_type == "Initial":
            initialDataType = str(params['initialDataType'])

            size = params["xScale"]
            resolution = int(params["resolution"])
            self.calc.resolution = resolution
            self.calc.size = size

            logging.info("Size: {} x {} Mpc^2, resolution: {} x {}".format(size, size, resolution, resolution))

            SIlimit = params['SILimit']

            if SIlimit == "None":
                SIlimit = None

            sigma = float(params['sigma'])

            SI_ns = params['n_s']
            if initialDataType == "SI":
                A_s = 2.214 * 10**(-9)
            else:
                A_s = 1

            redshift = generate_redshifts(params["redshift"])
            self.calc.redshift = redshift

            self.write_message(
                json.dumps({
                    'type': 'redshift',
                    'redshift': redshift.tolist()
                }))

            logging.info("Submitting initial state generation to ThreadPoolExecutor")
            yield pool.submit(self.set_initial_condition, sigma, initialDataType,
                              SIlimit, SI_ns, A_s)
            self.send_initial_state()
            self.write_message(json.dumps({'type': 'success', 'sort': 'Initial'}))

        elif param_type == "Cosmo":
            logging.info("Received cosmological parameters")
            cosmological_parameters = params
            logging.info("Submitting calculation to ThreadPoolExecutor")
            messages = yield pool.submit(self.set_cosmological_parameters, cosmological_parameters)
            for message in messages:
                self.write_message(json.dumps(message))
        elif param_type == "Start":
            logging.info("Starting propagation...")
            try:
                for redindex, z in enumerate(self.calc.redshift):
                    self.send_frame(redindex)
                self.write_message(json.dumps({'type': 'success', 'sort': 'Data'}))
            except Exception as e:
                logging.exception(e)
                self.send_exception(e)

    def send_frame(self, redindex):
        # `extrema`: (minimum, maximum) of (real space) data
        Valuenew, FValuenew, extrema = self.calc.getData(redindex)
        logging.info("Sending data for redshift = {}".format(self.calc.redshift[redindex]))

        # Create data to be displayed in transfer function window
        TransferData, _ = self.calc.getTransferData(redindex)

        self.write_message(json.dumps({'type': 'extrema', 'extrema': extrema}))
        progress = float(redindex) / len(self.calc.redshift)

        real = {quantity: base64.b64encode(data.astype(np.float32)) for quantity, data in Valuenew.iteritems()}
        transfer = {quantity: base64.b64encode(data.astype(np.float32)) for quantity, data in TransferData.iteritems()}
        self.write_message(
            json.dumps({
                'type': 'data',
                'progress': progress,
                'real': real,
                'fourier': [],
                'transfer': transfer,
            }))

    def send_initial_state(self):
        Value, FValue, extrema = self.calc.getInitialData()
        TransferData = np.ones(config.TRANSFER_FUNCTION_CLIENT_SAMPLES)
        krange = np.zeros(config.TRANSFER_FUNCTION_CLIENT_SAMPLES)
        logging.info("Sending initial data to client.")
        self.write_message({
            "type": "resolution",
            "value": self.calc.resolution
            })
        extremastring = json.dumps({'type': 'extrema', 'extrema': extrema})
        datastring = json.dumps({
            'type': 'data',
            'real': base64.b64encode(Value.astype(np.float32)),
            'fourier': [],
            'transfer': base64.b64encode(TransferData.astype(np.float32)),
            'k': krange.tolist()
            })
        self.write_message(extremastring)
        self.write_message(datastring)


    def set_initial_condition(self, sigma, initialDataType, SIlimit, SI_ns, A_s):
        try:
            self.calc.setInitialConditions(
                sigma=sigma,
                initialDataType=initialDataType,
                SIlimit=SIlimit,
                SI_ns=SI_ns,
                A=A_s
                )
        except Exception as e:
            logging.exception(e)
            self.send_exception(e)

    def send_exception(self, e):
        self.write_message(json.dumps({'type': 'exception', 'exception': traceback.format_exc()}))

    def set_cosmological_parameters(self, cosmologicalParameters):
        try:
            messages = []
            logging.info("Starting calculation...")
            self.calc.setCosmologialParameters(cosmologicalParameters=cosmologicalParameters)
            logging.info("Finished calculation!")

            messages.append({'type': 'success', 'sort': 'Cosmo'})
            messages.append({
                'type': 'Cl',
                'l': self.calc.tCl.l.tolist(),
                'tCl': self.calc.tCl.tCl.tolist()
                })
            messages.append({
                'type': 'mPk',
                'kh': self.calc.mPk.kh.tolist(),
                'Pkh': self.calc.mPk.Pkh.tolist()
                })

            z_of_decoupling = self.calc.z_dec
            frame_of_decoupling = np.argmin(np.abs(z_of_decoupling - self.calc.redshift))
            if self.calc.redshift[frame_of_decoupling] > z_of_decoupling:
                frame_of_decoupling -= 1
            messages.append({
                'type': 'decoupling',
                'frame': frame_of_decoupling,
                'z': z_of_decoupling})
        except Exception as e:
            logging.exception(e)
            self.send_exception(e)
        else:
            return messages


def main():
    logging.getLogger().setLevel(logging.DEBUG)

    application = tornado.web.Application(
        [
            (r"/", SimulationHandler),
            (r"/datasocket", DataConnection),
        ],
        template_path=os.path.join(os.path.dirname(__file__), "templates"),
        static_path=os.path.join(os.path.dirname(__file__), "static"),
        debug=True,
    )

    PORT = config.PORT if len(sys.argv) == 1 else int(sys.argv[1])
    application.listen(PORT)
    logging.info("Application launched on http://localhost:{}".format(PORT))
    IOLoop.instance().current().start()


if __name__ == '__main__':
    main()
