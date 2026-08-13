import logging
import vcf
from enum import Enum
_log = logging.getLogger(__name__)


class NoInputFileCoordinateSystemSelectedException(Exception):
    pass

class variant_coordinate_system(Enum):
    """For differences, see: https://www.biostars.org/p/84686/"""
    zero_based = 1
    one_based = 2

def tabix_query(filename, chrom, start, end, inputfile_variant_coordinate_system, encoding='utf-8'):
    """Call tabix and generate an array of strings for each line it returns."""
    vcf_reader = vcf.Reader(filename=filename, encoding=encoding)

    if inputfile_variant_coordinate_system == variant_coordinate_system.one_based:
        records = vcf_reader.fetch(chrom, start - 1, end)
    elif inputfile_variant_coordinate_system == variant_coordinate_system.zero_based:
        records = vcf_reader.fetch(chrom, start, end)
    else:
        raise NoInputFileCoordinateSystemSelectedException("No valid input file variant coordinate system was selected to adjust for when querying variant file")

    while True:
        try:
            record = next(records)
        except StopIteration:
            return
        except ValueError as e:
            # a single malformed record must not abort the remainder of the region
            _log.error("Error while parsing record from '" + filename + "' in region 'chr" + str(chrom) + ":" + str(start) + "-" + str(end) + "' :" + str(e))
            continue
        yield record