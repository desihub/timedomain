import os
import datetime

pi = 3.141592653589793


def radec_str2geojson(ra_str, dec_str):

    # hms -> ::, dms -> ::
    if isinstance(ra_str, str) and isinstance(dec_str, str):
        if ("h" in ra_str) and ("m" in ra_str) and ("s" in ra_str):
            ra_str = ra_str[:-1]  # strip 's' at the end
            for char in ("h", "m"):
                ra_str = ra_str.replace(char, ":")
        if ("d" in dec_str) and ("m" in dec_str) and ("s" in dec_str):
            dec_str = dec_str[:-1]  # strip 's' at the end
            for char in ("d", "m"):
                dec_str = dec_str.replace(char, ":")

        if (":" in ra_str) and (":" in dec_str):
            ra, dec = radec_str2rad(ra_str, dec_str)
            # convert to geojson-friendly degrees:
            ra = ra * 180.0 / pi - 180.0
            dec = dec * 180.0 / pi
        else:
            raise Exception("unrecognized string ra/dec format.")
    else:
        # already in degrees?
        ra = float(ra_str)
        # geojson-friendly ra:
        ra -= 180.0
        dec = float(dec_str)

    return ra, dec


def radec_str2rad(_ra_str, _dec_str):
    """
    :param _ra_str: 'H:M:S'
    :param _dec_str: 'D:M:S'
    :return: ra, dec in rad
    """
    # convert to rad:
    _ra = list(map(float, _ra_str.split(":")))
    _ra = (_ra[0] + _ra[1] / 60.0 + _ra[2] / 3600.0) * pi / 12.0
    _dec = list(map(float, _dec_str.split(":")))
    _sign = -1 if _dec_str.strip()[0] == "-" else 1
    _dec = (
        _sign
        * (abs(_dec[0]) + abs(_dec[1]) / 60.0 + abs(_dec[2]) / 3600.0)
        * pi
        / 180.0
    )

    return _ra, _dec


def load_config(path="/app", config_file="config.yaml"):
    """
    Load config and secrets
    """
    with open(os.path.join(path, config_file)) as cyaml:
        config = yaml.load(cyaml, Loader=yaml.FullLoader)

    return config


def log(message):
    print(f"{time_stamp()}: {message}")


def time_stamp():
    """
    :return: UTC time -> string
    """
    return datetime.datetime.utcnow().strftime("%Y%m%d_%H:%M:%S")


class Mongo:
    def __init__(
        self,
        host: str = "127.0.0.1",
        port: str = "27017",
        username: str = None,
        password: str = None,
        db: str = None,
        verbose=0,
        **kwargs,
    ):
        self.host = host
        self.port = port
        self.username = username
        self.password = password

        self.client = pymongo.MongoClient(host=self.host, port=self.port)
        self.db = self.client[db]
        # authenticate
        self.db.authenticate(self.username, self.password)

        self.verbose = verbose

    def insert_one(
        self, collection: str, document: dict, transaction: bool = False, **kwargs
    ):
        # note to future me: single-document operations in MongoDB are atomic
        # turn on transactions only if running a replica set
        try:
            if transaction:
                with self.client.start_session() as session:
                    with session.start_transaction():
                        self.db[collection].insert_one(document, session=session)
            else:
                self.db[collection].insert_one(document)
        except Exception as e:
            if self.verbose:
                print(
                    time_stamp(),
                    f"Error inserting document into collection {collection}: {str(e)}",
                )
                traceback.print_exc()

    def insert_many(
        self, collection: str, documents: list, transaction: bool = False, **kwargs
    ):
        ordered = kwargs.get("ordered", False)
        try:
            if transaction:
                with self.client.start_session() as session:
                    with session.start_transaction():
                        self.db[collection].insert_many(
                            documents, ordered=ordered, session=session
                        )
            else:
                self.db[collection].insert_many(documents, ordered=ordered)
        except BulkWriteError as bwe:
            if self.verbose:
                print(
                    time_stamp(),
                    f"Error inserting documents into collection {collection}: {str(bwe.details)}",
                )
                traceback.print_exc()
        except Exception as e:
            if self.verbose:
                print(
                    time_stamp(),
                    f"Error inserting documents into collection {collection}: {str(e)}",
                )
                traceback.print_exc()

    def update_one(
        self,
        collection: str,
        filt: dict,
        update: dict,
        transaction: bool = False,
        **kwargs,
    ):
        upsert = kwargs.get("upsert", True)

        try:
            if transaction:
                with self.client.start_session() as session:
                    with session.start_transaction():
                        self.db[collection].update_one(
                            filter=filt,
                            update=update,
                            upsert=upsert,
                            session=session,
                        )
            else:
                self.db[collection].update_one(
                    filter=filt, update=update, upsert=upsert
                )
        except Exception as e:
            if self.verbose:
                print(
                    time_stamp(),
                    f"Error inserting document into collection {collection}: {str(e)}",
                )
                traceback.print_exc()



