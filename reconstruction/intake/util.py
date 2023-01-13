import csv

def readClassificationCsv(csvFile):
	classifications = {}
	reader = csv.reader(csvFile)
	for row in reader:
		item, classification = row
		if not classification in classifications:
			classifications[classification] = []
		classifications[classification].append(item)
	return classifications

