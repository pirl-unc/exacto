from pkg_resources import get_distribution
from .variant_annotation import VariantAnnotation
from .variant_call import VariantCall
from .variant_filter import VariantFilter
from .variant import Variant
from .variants_list import VariantsList


__version__ = get_distribution('Exacto').version
__all__ = [
    'VariantAnnotation',
    'VariantCall',
    'VariantFilter',
    'Variant',
    'VariantsList'
]