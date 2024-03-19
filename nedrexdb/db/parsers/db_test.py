import pymongo
from configparser import ConfigParser



# Holen Sie sich die MongoDB-Konfiguration für die Entwicklungsdatenbank
mongo_port = 27019
mongo_db_name = "licensed_nedrex_dev"

# Verbinde dich mit der MongoDB-Datenbank
client = pymongo.MongoClient('localhost', mongo_port)
db = client[mongo_db_name]

# Beispielabfrage
collection = db['GeneAssociatedWithDisorder']
num_documents = collection.count_documents({})

print("Anzahl der Dokumente:", num_documents)
