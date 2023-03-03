from ..ConfigParser import *

from numbers import Number

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
	TypedConfigItem("networkTreatments", list),
	TypedConfigItem("correlationMethod", str, default="spearman"),
	TypedConfigItem("correlationCombinePValuesMethod", str, default="fisher"),
	TypedConfigItem("correlationCorrectPValuesMethod", str, default="fdr"),
	TypeOrSubconfigConfigItem("correlationPValueThresholds", Number, Config([
		MultiTypeConfigItem("individual", [Number, dict], default=0.2),
		MultiTypeConfigItem("combined", [Number, dict], default=0.05),
		MultiTypeConfigItem("corrected", [Number, dict], default=0.1)
	]), default=UseSubDefaults()),
	TypedConfigItem("correlationFilterMethod", str, default="allsamesign"),
	TypedConfigItem("correlationFilterPercentAgreementThreshold", Number, default=0.75),
	TypedConfigItem("correctCorrelationPValuesAfterConsistencyFiltering", bool, default=False)
])

