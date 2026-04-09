from importlib.metadata import PackageNotFoundError, version

from .codon_usage import CodonUsageTable
from .sequence_analyzer import SequenceAnalyzer
from .trna_weights import TRNAWeightTable

try:
	__version__ = version("seqmetrics")
except PackageNotFoundError:
	__version__ = "0.1.0"

__all__ = [
	"CodonUsageTable",
	"SequenceAnalyzer",
	"TRNAWeightTable",
]
