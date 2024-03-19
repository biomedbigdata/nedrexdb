import pymongo

# Verbindung zur MongoDB-Datenbank herstellen
client = pymongo.MongoClient('localhost', 27020)
db = client['licensed_nedrex_live']  # Datenbank auswählen

# Ein zufälliges Element aus der Sammlung abrufen
random_document = db['gene_associated_with_disorder'].find_one()

# Ergebnis ausgeben
print(random_document)
print(random_document)
