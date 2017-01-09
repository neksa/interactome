import csv
from collections import namedtuple, defaultdict
from itertools import islice
# from werkzeug.contrib.cache import SimpleCache
from sqlalchemy import create_engine, MetaData, Table


engine = create_engine('mysql+pymysql://gonceare:7lkj2Khsp42%l@localhost/interactomes', convert_unicode=True, echo=True)
metadata = MetaData(bind=engine)

species = Table('species', metadata, autoload=True)
pairs = Table('pairs', metadata, autoload=True)
proteins = Table('proteins', metadata, autoload=True)
idmapping = Table('idmapping', metadata, autoload=True)


def get_data():
    return []


if __name__ == '__main__':
    pass
