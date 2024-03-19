import pymongo
from configparser import ConfigParser



# Holen Sie sich die MongoDB-Konfiguration für die Entwicklungsdatenbank
mongo_port = 27020
mongo_db_name = "licensed_nedrex_live"

# Verbinde dich mit der MongoDB-Datenbank
client = pymongo.MongoClient('localhost', mongo_port)
db = client[mongo_db_name]

# Beispielabfrage
collections = db.list_collection_names()

# Ausgabe der Collections
for collection in collections:
    print(collection)
