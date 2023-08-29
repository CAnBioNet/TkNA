from ..ConfigParser import *

from numbers import Number

def validateConfig(config):
	metatreatments = config["metatreatments"]
	metatreatmentsToCombine = config["metatreatmentsToCombine"]
	metatreatmentsForDirectionFiltering = config["metatreatmentsForDirectionFiltering"]
	if metatreatments is not None:
		if config["networkTreatment"] is not None:
			raise ValidationFailedError("\"networkTreatment\" and \"metatreatments\" cannot both be defined. Either define \"networkTreatment\" to construct the network from a single treatment or define \"metatreatments\" to construct it from pooled treatments.")

		# Convert experiment/treatments to actual tuples (not supported by JSON)
		for metatreatmentName in metatreatments:
			metatreatments[metatreatmentName] = [tuple(experimentAndTreatment) for experimentAndTreatment in metatreatments[metatreatmentName]]
		config["metatreatments"] = metatreatments

		# Ensure metatreatments do not overlap
		experimentAndTreatmentsSeen = []
		for metatreatment in metatreatments.values():
			for experimentAndTreatment in metatreatment:
				if experimentAndTreatment in experimentAndTreatmentsSeen:
					raise ValidationFailedError(f"Metatreatments cannot overlap. Offending experiment/treatment: {experimentAndTreatment}")
				experimentAndTreatmentsSeen.append(experimentAndTreatment)

	elif config["networkTreatment"] is not None:
		if metatreatmentsToCombine is not None:
			raise ValidationFailedError("\"networkTreatment\" and \"metatreatmentsToCombine\" cannot both be defined. Either define \"networkTreatment\" only to construct the network from a single treatment or define \"metatreatments\" to construct it from pooled treatments.")
		elif metatreatmentsForDirectionFiltering is not None:
			raise ValidationFailedError("\"networkTreatment\" and \"metatreatmentsForDirectionFiltering\" cannot both be defined. Either define \"networkTreatment\" only to construct the network from a single treatment or define \"metatreatments\" to construct it from pooled treatments.")
	else:
		raise ValidationFailedError("At least one of \"networkTreatment\" and \"metatreatments\" must be defined. Either define \"networkTreatment\" to construct the network from a single treatment or define \"metatreatments\" to construct it from pooled treatments.")

	return config

# TODO: Can put descriptions here
singleCellConfigSpec = Config([
	TypedConfigItem("combineCellsMethod", str, "mean"),
	TypedConfigItem("differenceCombinePValuesMethod", str, default="fisher"),
	TypeOrSubconfigConfigItem("differencePValueThresholds", Number, Config([
		MultiTypeConfigItem("individual", [Number, dict], default=0.2),
		MultiTypeConfigItem("combined", [Number, dict], default=0.05),
		MultiTypeConfigItem("corrected", [Number, dict], default=0.1)
	]), default=UseSubDefaults()),
	TypedConfigItem("foldChangeFilterMethod", str, default="allsamesign"),
	TypedConfigItem("foldChangeFilterPercentAgreementThreshold", Number, default=0.75),
	TypedConfigItem("networkTreatment", str, default=None),
	TypedConfigItem("correlationMethod", str, default="spearman"),
	TypedConfigItem("correlationCombinePValuesMethod", str, default="fisher"),
	TypedConfigItem("correlationCorrectPValuesMethod", str, default="fdr"),
	TypeOrSubconfigConfigItem("correlationPValueThresholds", Number, Config([
		MultiTypeConfigItem("individual", [Number, dict], default=0.2),
		MultiTypeConfigItem("combined", [Number, dict], default=0.05),
		MultiTypeConfigItem("corrected", [Number, dict], default=0.1)
	]), default=UseSubDefaults()),
	MultiTypeConfigItem("correlationCoefficientThresholds", [Number, dict], default=None),
	TypedConfigItem("correlationFilterMethod", str, default="allsamesign"),
	TypedConfigItem("correlationFilterPercentAgreementThreshold", Number, default=0.75),
	TypedConfigItem("correctCorrelationPValuesAfterConsistencyFiltering", bool, default=False),
	TypedConfigItem("metatreatments", dict, default=None),
	TypedConfigItem("metatreatmentsToCombine", list, default=None),
	TypedConfigItem("metatreatmentsForDirectionFiltering", list, default=None),
	TypedConfigItem("noPUC", bool, default=False)
], validateConfig)

