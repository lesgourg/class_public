import pickle
import os
import logging
import uuid

class Database:
    def __init__(self, directory, db_file="database.dat"):
        self.directory = directory
        self.db_file = db_file

        if not os.path.isdir(directory):
            raise ValueError("'{}' is not a directory!".format(directory))

        self.db_path = os.path.join(directory, db_file)
        if not os.path.exists(self.db_path):
            logging.info("No database found; Creating one at {}.".format(self.db_path))
            with open(self.db_path, "w") as f:
                pickle.dump(dict(), f)

        self.db = self.__read_database()

    def __read_database(self):
        with open(self.db_path) as f:
            return pickle.load(f)

    def __write_database(self):
        with open(self.db_path, "w") as f:
            pickle.dump(self.db, f)

    def __create_file(self, data):
        filename = str(uuid.uuid4())
        with open(os.path.join(self.directory, filename), "w") as f:
            pickle.dump(data, f)
        return filename

    def __get_frozen_key(self, key):
        return frozenset(key.items())

    def __getitem__(self, key):
        frozen_key = self.__get_frozen_key(key)
        if frozen_key in self.db:
            filename = self.db[frozen_key]
            with open(os.path.join(self.directory, filename)) as f:
                return pickle.load(f)
        else:
            raise KeyError("No data for key: {}".format(key))

    def __setitem__(self, key, data):
        frozen_key = self.__get_frozen_key(key)
        self.db[frozen_key] = self.__create_file(data)
        self.__write_database()

    def __contains__(self, key):
        """
        Return whether `self` contains a record
        for the given `key`.
        """
        return self.__get_frozen_key(key) in self.db