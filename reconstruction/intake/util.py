import chardet
import csv

def readAndDecodeFile(filePath):
	with open(filePath, "rb") as targetFile:
		fileBytes = targetFile.read()
		fileEncoding = chardet.detect(fileBytes)["encoding"]
		return fileBytes.decode(fileEncoding)

def readClassificationCsv(csvString):
	classifications = {}
	reader = csv.reader(csvString.splitlines())
	for row in reader:
		item, classification = row
		item = item.strip()
		classification = classification.strip()
		if not classification in classifications:
			classifications[classification] = []
		classifications[classification].append(item)
	return classifications

